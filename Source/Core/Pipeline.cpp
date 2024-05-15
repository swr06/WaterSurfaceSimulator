#include "Pipeline.h"

#include "Utils/Random.h"

#include "Object.h"

#include "FpsCamera.h"
#include "Player.h"

namespace Simulation {

	float CurrentTime = glfwGetTime();
	float Frametime = 0.0f;
	float DeltaTime = 0.0f;

	Player MainPlayer;
	FPSCamera& Camera = MainPlayer.Camera;

	std::vector<Object> Objects;

	class RayTracerApp : public Simulation::Application
	{
	public:

		bool vsync;

		RayTracerApp()
		{
			m_Width = 800;
			m_Height = 600;
		}

		void OnUserCreate(double ts) override
		{

		}

		void OnUserUpdate(double ts) override
		{
			glfwSwapInterval((int)vsync);

			GLFWwindow* window = GetWindow();

		}

		void OnImguiRender(double ts) override
		{
			ImGuiIO& io = ImGui::GetIO();
			if (ImGui::Begin("Debug/Edit Mode")) {

				ImGui::Text("Camera Position : %f,  %f,  %f", Camera.GetPosition().x, Camera.GetPosition().y, Camera.GetPosition().z);
				ImGui::Text("Camera Front : %f,  %f,  %f", Camera.GetFront().x, Camera.GetFront().y, Camera.GetFront().z);
				ImGui::Text("Time : %f s", glfwGetTime());
			
			} ImGui::End();
		}

		void OnEvent(Simulation::Event e) override
		{
			ImGuiIO& io = ImGui::GetIO();

			if (e.type == Simulation::EventTypes::MousePress && !ImGui::GetIO().WantCaptureMouse && GetCurrentFrame() > 32)
			{

			}

			if (e.type == Simulation::EventTypes::MouseMove && GetCursorLocked())
			{
				Camera.UpdateOnMouseMovement(e.mx, e.my);
			}


			if (e.type == Simulation::EventTypes::MouseScroll && !ImGui::GetIO().WantCaptureMouse)
			{
				float Sign = e.msy < 0.0f ? 1.0f : -1.0f;
				Camera.SetFov(Camera.GetFov() + 2.0f * Sign);
				Camera.SetFov(glm::clamp(Camera.GetFov(), 1.0f, 89.0f));
			}

			if (e.type == Simulation::EventTypes::WindowResize)
			{
				Camera.SetAspect((float)glm::max(e.wx, 1) / (float)glm::max(e.wy, 1));
			}

			if (e.type == Simulation::EventTypes::KeyPress && e.key == GLFW_KEY_ESCAPE) {
				exit(0);
			}

			if (e.type == Simulation::EventTypes::KeyPress && e.key == GLFW_KEY_F1)
			{
				this->SetCursorLocked(!this->GetCursorLocked());
			}

			if (e.type == Simulation::EventTypes::KeyPress && e.key == GLFW_KEY_F2 && this->GetCurrentFrame() > 5)
			{
				Simulation::ShaderManager::RecompileShaders();
			}

			if (e.type == Simulation::EventTypes::KeyPress && e.key == GLFW_KEY_F3 && this->GetCurrentFrame() > 5)
			{
				Simulation::ShaderManager::ForceRecompileShaders();
			}

			if (e.type == Simulation::EventTypes::KeyPress && e.key == GLFW_KEY_V && this->GetCurrentFrame() > 5)
			{
				vsync = !vsync;
			}



		}


	};

	void Pipeline::StartPipeline()
	{
		// Application
		RayTracerApp app;
		app.Initialize();
		app.SetCursorLocked(false);

		// Create VBO and VAO for drawing the screen-sized quad.
		GLClasses::VertexBuffer ScreenQuadVBO;
		GLClasses::VertexArray ScreenQuadVAO;

		GLClasses::VertexBuffer WaterMeshVBO;
		GLClasses::IndexBuffer WaterMeshEBO;
		GLClasses::VertexArray WaterMeshVAO;

		GLClasses::Texture Heightmap;
		Heightmap.CreateTexture("Res/Heightmap.png", false, false, false, GL_TEXTURE_2D,
			GL_LINEAR, GL_LINEAR,
			GL_REPEAT, GL_REPEAT, false);

		// Setup screensized quad for rendering
		{
			unsigned long long CurrentFrame = 0;
			float QuadVertices_NDC[] =
			{
				-1.0f,  1.0f,  0.0f, 1.0f, -1.0f, -1.0f,  0.0f, 0.0f,
				 1.0f, -1.0f,  1.0f, 0.0f, -1.0f,  1.0f,  0.0f, 1.0f,
				 1.0f, -1.0f,  1.0f, 0.0f,  1.0f,  1.0f,  1.0f, 1.0f
			};

			ScreenQuadVAO.Bind();
			ScreenQuadVBO.Bind();
			ScreenQuadVBO.BufferData(sizeof(QuadVertices_NDC), QuadVertices_NDC, GL_STATIC_DRAW);
			ScreenQuadVBO.VertexAttribPointer(0, 2, GL_FLOAT, 0, 4 * sizeof(GLfloat), 0);
			ScreenQuadVBO.VertexAttribPointer(1, 2, GL_FLOAT, 0, 4 * sizeof(GLfloat), (void*)(2 * sizeof(GLfloat)));
			ScreenQuadVAO.Unbind();
		}

		int Resolution = 128;

		{
			std::vector<float> BufferData;
			std::vector<unsigned int> Indices;


			float Range = 128.0f;

			for (int x = -Resolution; x <= Resolution; x++) {

				for (int y = -Resolution; y <= Resolution; y++) {

					float Cx = (Range / (float)(Resolution)) * float(x);
					float Cy = (Range / (float)(Resolution)) * float(y);

					BufferData.push_back(Cx);
					BufferData.push_back(Cy);
				}
			}

			int Height = Resolution * 2;

			for (int x = 0; x <= Height; x++) {
				for (int y = 0; y < Height; y++) {
					int i1 = (x * Height) + (Height + 2 + x + y) - 1;
					int i2 = (x * Height) + y + x + 1 - 1;
					int i3 = (x * Height) + y + x + 2 - 1;
					int i4 = (x * Height) + y + x + 2 - 1;
					int i5 = (x * Height) + (Height + 2 + y + x) + 1 - 1;
					int i6 = i1;

					Indices.push_back(i1);
					Indices.push_back(i2);
					Indices.push_back(i3);
					Indices.push_back(i4);
					Indices.push_back(i5);
					Indices.push_back(i6);
				}
			}

			WaterMeshVAO.Bind();
			WaterMeshVBO.Bind();
			WaterMeshEBO.Bind();
			WaterMeshVBO.BufferData(BufferData.size() * sizeof(float), BufferData.data(), GL_STATIC_DRAW);
			WaterMeshVBO.VertexAttribPointer(0, 2, GL_FLOAT, 0, 2 * sizeof(GLfloat), 0);
			WaterMeshEBO.BufferData(Indices.size() * sizeof(unsigned int), Indices.data(), GL_STATIC_DRAW);
			WaterMeshVAO.Unbind();
		}

		// Create Shaders 
		ShaderManager::CreateShaders();

		// Shaders
		GLClasses::Shader& BlitShader = ShaderManager::GetShader("BLIT");
		GLClasses::Shader& BasicRender = ShaderManager::GetShader("BASICRENDER");


		
		glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
		glClear(GL_COLOR_BUFFER_BIT);

		while (!glfwWindowShouldClose(app.GetWindow())) {

			glDisable(GL_DEPTH_TEST);
			glDisable(GL_CULL_FACE);

			app.OnUpdate();

			// Player
			MainPlayer.OnUpdate(app.GetWindow(), DeltaTime, 0.4f, app.GetCurrentFrame());

			glBindFramebuffer(GL_FRAMEBUFFER,0);
			glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
			glClear(GL_COLOR_BUFFER_BIT);
			glPolygonMode(GL_FRONT, GL_LINE);
			glPolygonMode(GL_BACK, GL_LINE);

			BasicRender.Use();
			BasicRender.SetMatrix4("u_ViewProj", Camera.GetViewProjection());
			BasicRender.SetMatrix4("u_ViewProjRot", Camera.GetViewProjection() * glm::rotate(glm::mat4(1.0f), 1.570f, glm::vec3(1.0f, 0.0f, 0.0f)));

			WaterMeshVAO.Bind();
			glDrawElements(GL_TRIANGLES, Resolution * Resolution * 4 * 6, GL_UNSIGNED_INT, 0);
			WaterMeshVAO.Unbind();

			glUseProgram(0);

			glFinish();
			app.FinishFrame();

			CurrentTime = glfwGetTime();
			DeltaTime = CurrentTime - Frametime;
			Frametime = glfwGetTime();

			GLClasses::DisplayFrameRate(app.GetWindow(), "Simulation ");
		}
	}
}