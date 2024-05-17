#include "Pipeline.h"

#include "Utils/Random.h"

#include "Object.h"

#include "FpsCamera.h"
#include "Player.h"

#define MAX_SPHERES 16

namespace Simulation {

	int Resolution = 192;
	float Range = 4.0f;
	float Width = Range / Resolution;

	int Substeps = 3;

	float c = 10.0f;
	float s = 1.0f;
	float kProportionality = (c * c) / (s * s);
	bool DoSim = false;
	float DampingCoeff = 0.4f;

	float* Heightmap;
	float* ObjectHeights;
	float* WaterVelocities;
	float* WaterAccelerations;
	Random RandomGen;

	int To1DIdx(int x, int y) {
		return (y * Resolution) + x;
	}

	struct Sphere {
		glm::vec3 Position;
		glm::vec3 Velocity;
		glm::vec3 Acceleration;
		float Mass;
		float Density;
		float Radius;
	};

	struct RenderSphere {
		glm::vec4 PositionRadius;
	};

	std::vector<Sphere> Spheres;

	float CurrentTime = glfwGetTime();
	float Frametime = 0.0f;
	float DeltaTime = 0.0f;
	bool WireFrame = false;

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
			static float r = 0.5f;

			ImGuiIO& io = ImGui::GetIO();
			if (ImGui::Begin("Debug/Edit Mode")) {

				ImGui::Text("Camera Position : %f,  %f,  %f", Camera.GetPosition().x, Camera.GetPosition().y, Camera.GetPosition().z);
				ImGui::Text("Camera Front : %f,  %f,  %f", Camera.GetFront().x, Camera.GetFront().y, Camera.GetFront().z);
				ImGui::Text("Time : %f s", glfwGetTime());

				ImGui::NewLine();
				ImGui::Checkbox("Do Sim", &DoSim);
				ImGui::SliderInt("Substeps", &Substeps, 1, 100);
				ImGui::SliderFloat("c", &c, 0.0f, 100.0f);
				ImGui::SliderFloat("Damping Coeff", &DampingCoeff, 0.01f, 4.0f);
				ImGui::SliderFloat("s", &s, 0.1f, 100.0f);

				if (ImGui::Button("Mod")) {
					Heightmap[int(RandomGen.Float() * Resolution * Resolution)] = 1.5;
					Heightmap[int(RandomGen.Float() * Resolution * Resolution)] = 0.5f;
				}
				
				if (ImGui::Button("Reset")) {
					
					for (int i = 0; i < Resolution * Resolution; i++) {
						Heightmap[i] = 1.f;
					}

					memset(WaterAccelerations, 0, Resolution * Resolution * sizeof(float));
					memset(WaterVelocities, 0, Resolution * Resolution * sizeof(float));
					memset(ObjectHeights, 0, Resolution * Resolution * sizeof(float));

				}



				ImGui::NewLine();

				ImGui::Checkbox("Wireframe", &WireFrame);

				ImGui::NewLine();
				ImGui::NewLine();

				if (Spheres.size() < MAX_SPHERES) {
					if (ImGui::Button("Place Sphere")) {
						Sphere s1 = { glm::vec3(Camera.GetPosition() + Camera.GetFront() * (r * 2.1f)), glm::vec3(0.0f), glm::vec3(0.0f), 1.0f, 1.0f, r };
						Spheres.push_back(s1);
					}

					ImGui::SliderFloat("Radii of Sphere", &r, 0.05f, 2.0f);
				}
				ImGui::NewLine();

				if (Spheres.size() > 0) {
					if (ImGui::Button("Delete Sphere")) {
						Spheres.pop_back();
					}
				}

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

	float SampleHeight(int x, int y) {
		return Heightmap[To1DIdx(x, y)];
	}

	glm::vec2 ConvertToWorldSpace(const glm::ivec2& Texel) {
		return ((glm::vec2(Texel) / float(Resolution)) * 2.0f - 1.0f) * Range;
	}

	float SampleHeightClamped(int x, int y) {
		x = glm::clamp(x, 0, Resolution - 1);
		y = glm::clamp(y, 0, Resolution - 1);
		return SampleHeight(x, y);
	}

	float SampleHeightClamped(glm::ivec2 t) {
		t.x = glm::clamp(t.x, 0, Resolution - 1);
		t.y = glm::clamp(t.y, 0, Resolution - 1);
		return SampleHeight(t.x, t.y);
	}

	void SimulateWaterAcceleration(float Dt) {

		glm::ivec2 NeighbourOffsets[4] = { glm::ivec2(-1,0), glm::ivec2(1,0), glm::ivec2(0, -1), glm::ivec2(0, 1) };

		for (int x = 0; x < Resolution; x++) {
			for (int y = 0; y < Resolution; y++) {

				float H = SampleHeightClamped(x, y);

				float Neighbours[4];

				float SigmaH = 0.;

				for (int i = 0; i < 4; i++) {
					Neighbours[i] = SampleHeightClamped(glm::ivec2(x, y) + NeighbourOffsets[i]);
					SigmaH += Neighbours[i];
				}

				float Acceleration = kProportionality * (SigmaH - 4.0f * H);

				WaterAccelerations[To1DIdx(x, y)] = Acceleration;

			}
		}

	}

	void SimulateWaterVelocities(float Dt) {

		for (int x = 0; x < Resolution; x++) {
			for (int y = 0; y < Resolution; y++) {

				float Acceleration = WaterAccelerations[To1DIdx(x, y)];

				WaterVelocities[To1DIdx(x, y)] += Acceleration * Dt;
				WaterVelocities[To1DIdx(x, y)] *= glm::clamp(1.0f - (DampingCoeff * Dt), 0.0f, 1.0f);
			}
		}

	}

	void SimulateWaterHeights(float Dt) {

		for (int x = 0; x < Resolution; x++) {
			for (int y = 0; y < Resolution; y++) {

				float Velocity = WaterVelocities[To1DIdx(x, y)];
				Heightmap[To1DIdx(x, y)] += Velocity * Dt;
			}
		}

	}

	void SimulateWater(float Dt) {



		float Ddt = Dt / float(Substeps);

		for (int i = 0; i < Substeps; i++) {

			SimulateWaterAcceleration(Ddt);
			SimulateWaterVelocities(Ddt);
			SimulateWaterHeights(Ddt);
		}

	}

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

		RandomGen.Float();

		//GLClasses::Texture Heightmap;
		//Heightmap.CreateTexture("Res/Heightmap.png", false, false, false, GL_TEXTURE_2D,
		//	GL_LINEAR, GL_LINEAR,
		//	GL_REPEAT, GL_REPEAT, false);

		// Setup screensized quad for rendering
		{
			float QuadVertices_NDC[] =
			{
				-1.0f,  1.0f,  0.0f, 1.0f, -1.0f, -1.0f,  0.0f, 0.0f,
				 1.0f, -1.0f,  1.0f, 0.0f, -1.0f,  1.0f,  0.0f, 1.0f,
				 1.0f, -1.0f,  1.0f, 0.0f,  1.0f,  1.0f,  1.0f, 1.0f
			};

			ScreenQuadVBO.Bind();
			ScreenQuadVBO.BufferData(sizeof(QuadVertices_NDC), QuadVertices_NDC, GL_STATIC_DRAW);
			
			ScreenQuadVAO.Bind();
			ScreenQuadVBO.Bind();
			ScreenQuadVBO.VertexAttribPointer(0, 2, GL_FLOAT, 0, 4 * sizeof(GLfloat), 0);
			ScreenQuadVBO.VertexAttribPointer(1, 2, GL_FLOAT, 0, 4 * sizeof(GLfloat), (void*)(2 * sizeof(GLfloat)));
			ScreenQuadVAO.Unbind();
		}


		{
			std::vector<float> BufferData;
			std::vector<unsigned int> Indices;

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
		GLClasses::Shader& RTSphere = ShaderManager::GetShader("SPHERE");

		GLClasses::Framebuffer GBuffer[2] = { GLClasses::Framebuffer(16, 16, {{GL_RGBA, GL_RGBA, GL_UNSIGNED_BYTE, false, false},  {GL_RGBA16F, GL_RGBA, GL_FLOAT, false, false}}, true, true), GLClasses::Framebuffer(16, 16, {{GL_RGBA, GL_RGBA, GL_UNSIGNED_BYTE, false, false}, {GL_RGBA16F, GL_RGBA, GL_FLOAT, false, false}}, true, true) };
		
		// Spheres
		GLuint SphereSSBO = 0;
		glGenBuffers(1, &SphereSSBO);
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, SphereSSBO);
		glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(RenderSphere) * MAX_SPHERES, (void*)0, GL_DYNAMIC_DRAW);
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

		// Create Heightmaps
		Heightmap = new float[Resolution * Resolution];
		ObjectHeights = new float[Resolution * Resolution];
		WaterVelocities = new float[Resolution * Resolution];
		WaterAccelerations = new float[Resolution * Resolution];
		memset(Heightmap, 0, Resolution* Resolution * sizeof(float));
		memset(WaterAccelerations, 0, Resolution* Resolution * sizeof(float));
		memset(WaterVelocities, 0, Resolution* Resolution * sizeof(float));
		memset(ObjectHeights, 0, Resolution* Resolution * sizeof(float));

		// GPU Data
		GLuint HeightmapSSBO = 0;
		glGenBuffers(1, &HeightmapSSBO);
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, HeightmapSSBO);
		glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * Resolution * Resolution, (void*)0, GL_DYNAMIC_DRAW);
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

		// For now.
		for (int x = 0; x < Resolution; x++) {
			for (int y = 0; y < Resolution; y++) {
				int i = To1DIdx(x, y);
				Heightmap[i] = 1.0f;
			}
		}

		// Make Spheres

		Sphere s1 = { glm::vec3(0.0f), glm::vec3(0.0f), glm::vec3(0.0f), 1.0f, 1.0f, 0.5f };
		Spheres.push_back(s1);

		glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
		glClear(GL_COLOR_BUFFER_BIT);

		while (!glfwWindowShouldClose(app.GetWindow())) {

			kProportionality = (c * c) / (s * s);

			glDisable(GL_CULL_FACE);

			app.OnUpdate();

			// FBO Update
			GBuffer[0].SetSize(app.GetWidth(), app.GetHeight());
			GBuffer[1].SetSize(app.GetWidth(), app.GetHeight());

			// Player
			MainPlayer.OnUpdate(app.GetWindow(), DeltaTime, 0.5f, app.GetCurrentFrame());

			// SIMULATE
			if (DoSim)
				SimulateWater(DeltaTime);

			///

			// Upadate spheres
			if (Spheres.size() > MAX_SPHERES) {
				throw "yo.";
			}

			std::vector<RenderSphere> Data;
			for (int i = 0; i < Spheres.size(); i++) {
				Data.push_back({ glm::vec4(Spheres[i].Position, Spheres[i].Radius) });
			}

			// Upload data
			glBindBuffer(GL_SHADER_STORAGE_BUFFER, SphereSSBO);
			glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(RenderSphere) * Data.size(), Data.data(), GL_DYNAMIC_DRAW);
			glBindBuffer(GL_SHADER_STORAGE_BUFFER, HeightmapSSBO);
			glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float)*Resolution*Resolution, Heightmap, GL_DYNAMIC_DRAW);

			// GBuffer
			GBuffer[0].Bind();
			glEnable(GL_DEPTH_TEST);
			glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			glPolygonMode(GL_FRONT, WireFrame ? GL_LINE : GL_FILL);
			glPolygonMode(GL_BACK, WireFrame ? GL_LINE : GL_FILL);

			BasicRender.Use();
			BasicRender.SetMatrix4("u_ViewProj", Camera.GetViewProjection());
			BasicRender.SetMatrix4("u_ModelMatrix", glm::rotate(glm::mat4(1.0f), 1.570f, glm::vec3(1.0f, 0.0f, 0.0f)));
			BasicRender.SetMatrix4("u_ViewProjRot", Camera.GetViewProjection() * glm::rotate(glm::mat4(1.0f), 1.570f, glm::vec3(1.0f, 0.0f, 0.0f)));
			BasicRender.SetFloat("u_Range", Range);
			BasicRender.SetFloat("u_RangeV", Range);
			BasicRender.SetInteger("u_Res", Resolution);
			BasicRender.SetInteger("u_ResV", Resolution);

			glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, HeightmapSSBO);

			WaterMeshVAO.Bind();
			glDrawElements(GL_TRIANGLES, Resolution * Resolution * 4 * 6, GL_UNSIGNED_INT, 0);
			WaterMeshVAO.Unbind();

			glDisable(GL_DEPTH_TEST);
			glDisable(GL_CULL_FACE);
			glUseProgram(0);
			glPolygonMode(GL_FRONT, GL_FILL);
			glPolygonMode(GL_BACK, GL_FILL);

			// Sphere

			GBuffer[1].Bind();
			RTSphere.Use();
			
			RTSphere.SetInteger("u_Texture", 0);
			RTSphere.SetInteger("u_Depth", 1);
			RTSphere.SetInteger("u_Spheres", Spheres.size());
			
			RTSphere.SetFloat("u_zNear", Camera.GetNearPlane());
			RTSphere.SetFloat("u_zFar", Camera.GetFarPlane());
			RTSphere.SetMatrix4("u_InverseProjection", glm::inverse(Camera.GetProjectionMatrix()));
			RTSphere.SetMatrix4("u_InverseView", glm::inverse(Camera.GetViewMatrix()));
			
			glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, SphereSSBO);
			
			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, GBuffer[0].GetTexture());
			
			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, GBuffer[0].GetDepthBuffer());
			
			ScreenQuadVAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			ScreenQuadVAO.Unbind();

			// Blit

			glBindFramebuffer(GL_FRAMEBUFFER, 0);

			glPolygonMode(GL_FRONT, GL_FILL);
			glPolygonMode(GL_BACK, GL_FILL);

			BlitShader.Use();

			BlitShader.SetInteger("u_Texture", 0);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, GBuffer[1].GetTexture());

			ScreenQuadVAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			ScreenQuadVAO.Unbind();

			glFinish();
			app.FinishFrame();

			CurrentTime = glfwGetTime();
			DeltaTime = CurrentTime - Frametime;
			Frametime = glfwGetTime();

			GLClasses::DisplayFrameRate(app.GetWindow(), "Simulation ");
		}
	}
}