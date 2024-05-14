#include "Pipeline.h"

#include "Utils/Random.h"

#include "Object.h"

namespace Simulation {

	float CurrentTime = glfwGetTime();
	float Frametime = 0.0f;
	float DeltaTime = 0.0f;

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
			}


			if (e.type == Simulation::EventTypes::MouseScroll && !ImGui::GetIO().WantCaptureMouse)
			{
			}

			if (e.type == Simulation::EventTypes::WindowResize)
			{
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

		// Create Shaders 
		ShaderManager::CreateShaders();

		// Shaders
		GLClasses::Shader& BlitShader = ShaderManager::GetShader("BLIT");
		float OrthographicRange = 400.0f;
		OrthographicCamera Orthographic(-400.0f, 400.0f, -400.0f, 400.0f);
		
		glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
		glClear(GL_COLOR_BUFFER_BIT);

		while (!glfwWindowShouldClose(app.GetWindow())) {

			glDisable(GL_DEPTH_TEST);
			glDisable(GL_CULL_FACE);

			app.OnUpdate();
			glBindFramebuffer(GL_FRAMEBUFFER,0);
			glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
			glClear(GL_COLOR_BUFFER_BIT);

			BlitShader.Use();
			ScreenQuadVAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			ScreenQuadVAO.Unbind();

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