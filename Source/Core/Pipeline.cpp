#include "Pipeline.h"

#include "Utils/Random.h"

#include "Object.h"

#include "FpsCamera.h"
#include "Player.h"

#include "GLClasses/CubeTextureMap.h"

#define MAX_SPHERES 24

// in meters
#define MAX_WAVE_HEIGHT 8 

namespace Simulation {

	// TODO ; WEight for heightmap contrib

	const glm::vec2 PoissonDisk[32] = {
		glm::vec2(-0.613392, 0.617481), glm::vec2(0.751946, 0.453352),
		glm::vec2(0.170019, -0.040254), glm::vec2(0.078707, -0.715323),
		glm::vec2(-0.299417, 0.791925), glm::vec2(-0.075838, -0.529344),
		glm::vec2(0.645680, 0.493210), glm::vec2(0.724479, -0.580798),
		glm::vec2(-0.651784, 0.717887), glm::vec2(0.222999, -0.215125),
		glm::vec2(0.421003, 0.027070), glm::vec2(-0.467574, -0.405438),
		glm::vec2(-0.817194, -0.271096), glm::vec2(-0.248268, -0.814753),
		glm::vec2(-0.705374, -0.668203), glm::vec2(0.354411, -0.887570),
		glm::vec2(0.977050, -0.108615), glm::vec2(0.175817, 0.382366),
		glm::vec2(0.063326, 0.142369), glm::vec2(0.487472, -0.063082),
		glm::vec2(0.203528, 0.214331), glm::vec2(-0.084078, 0.898312),
		glm::vec2(-0.667531, 0.326090), glm::vec2(0.488876, -0.783441),
		glm::vec2(-0.098422, -0.295755), glm::vec2(0.470016, 0.217933),
		glm::vec2(-0.885922, 0.215369), glm::vec2(-0.696890, -0.549791),
		glm::vec2(0.566637, 0.605213), glm::vec2(-0.149693, 0.605762),
		glm::vec2(0.039766, -0.396100), glm::vec2(0.034211, 0.979980)
	};


	// kg/m^3
	const float RhoWater = 998.2f;
	const float RhoRubber = 1522.0f;
	const float RhoWood = 700.0f; // Cherrywood? or something.

	float DebugVar = 0.0f;

	int Resolution = 256;

	float MouseRippleSize = 0.015f;

	float Range = 4.0f;
	float PoolRange = 4.0f;
	float PoolHeight = 2.0f;
	float ColumnWidth = (2.0f * Range) / float(Resolution);

	float AlphaO = 0.4f;

	int Substeps = 1;

	int SMOOTHING_ITR = 0;

	float CurrDens = RhoRubber;
	float SphereFrictionCoefficient = 0.025f;
	float SphereVel = 4.0f;
	static float SphereRad = 0.25f;
	float c = 12.0f;
	float s = 1.0f;
	float kProportionality = (c * c) / (s * s);
	bool DoSim = false;
	float DampingCoeff = 0.25f;

	bool RenderSpheres = true;
	bool RenderPool = true;

	float WaterBlueness = 1.75f;

	bool DoContainerCollisions = false;
	bool DoSSCollisions = true;
	
	bool PhysicsStep = false;

	int CheckerStep = 0; // Global variable that oscillates b/w 0 and 1

	float* Heightmap;
	float* ObjectHeights[2];
	float* WaterVelocities;
	float* WaterAccelerations;
	Random RandomGen;

	typedef glm::vec3 Force;

	inline int To1DIdx(int x, int y) {
		return (y * Resolution) + x;
	}

	inline int To1DIdxSafe(int x, int y) {
		return glm::clamp((y * Resolution) + x, 0, (Resolution * Resolution) - 1);
	}

	class Sphere {
	public: 
		glm::vec3 Position;
		glm::vec3 Velocity;
		float Density;
		float Radius;
		glm::vec3 NetAcceleration;
		glm::vec3 NetForce;

		inline float Mass() const noexcept {
			const float m = (4.0f / 3.0f) * 3.141592653;
			return m * Radius * Radius * Radius * Density;
		}

		inline float Volume() const noexcept {
			const float m = (4.0f / 3.0f) * 3.141592653;
			return m * Radius * Radius * Radius;
		}
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

	glm::vec3 RSI(glm::vec3 Origin, glm::vec3 Dir, float Radius)
	{
		using namespace glm;
		float VoV = dot(Dir, Dir);
		float Acc = VoV * Radius * Radius;

		// Solve quadratic 
		Acc += 2.0 * Origin.x * dot(glm::vec2(Origin.y, Origin.z), glm::vec2(Dir.y, Dir.z)) * Dir.x;
		Acc += 2.0 * Origin.y * Origin.z * Dir.y * Dir.z;
		Acc -= dot(Origin * Origin, vec3(dot(glm::vec2(Dir.y, Dir.z), glm::vec2(Dir.y, Dir.z)), dot(glm::vec2(Dir.x, Dir.z), glm::vec2(Dir.x, Dir.z)), dot(glm::vec2(Dir.x, Dir.y), glm::vec2(Dir.x, Dir.y))));

		// No intersect 
		if (Acc < 0.0)
		{
			return glm::vec3(-1.0f);
		}

		Acc = sqrt(Acc);

		float Dist1 = (Acc - dot(Origin, Dir)) / VoV;
		float Dist2 = -(Acc + dot(Origin, Dir)) / VoV;

		if (Dist1 >= 0.0 && Dist2 >= 0.0)
		{
			return glm::vec3(Dist1, Dist2, min(Dist1, Dist2));
		}
		else
		{
			return glm::vec3(Dist1, Dist2, max(Dist1, Dist2));
		}
	}

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

				ImGui::SliderFloat("DebugVar", &DebugVar, -1., 1.0f);
				ImGui::NewLine();

				ImGui::Text("Camera Position : %f,  %f,  %f", Camera.GetPosition().x, Camera.GetPosition().y, Camera.GetPosition().z);
				ImGui::Text("Camera Front : %f,  %f,  %f", Camera.GetFront().x, Camera.GetFront().y, Camera.GetFront().z);
				ImGui::Text("Time : %f s", glfwGetTime());
				ImGui::SliderFloat("Speed", &MainPlayer.Speed, 0.01f, 5.0f);

				if (ImGui::Button("Reset Speed")) {
					MainPlayer.Speed = 0.5f;
				}

				ImGui::NewLine();
				ImGui::SliderFloat("Water Blueness", &WaterBlueness, 0.01f, 6.0f);
				ImGui::NewLine();
				ImGui::SliderFloat("Range of water volume", &Range, 0.1f, 32.0f);
				ImGui::SliderFloat("Range of Pool", &PoolRange, 0.1f, Range);
				ImGui::SliderFloat("Height of Pool", &PoolHeight, 0.1f, 32.0f);
				if (ImGui::Button("Snap Pool")) {
					PoolRange = Range;
				}
				ImGui::NewLine();

				ImGui::Checkbox("Render Spheres", &RenderSpheres);
				ImGui::Checkbox("Render Pool", &RenderPool);

				ImGui::NewLine();

				ImGui::Checkbox("Do Sim", &DoSim);
				ImGui::Checkbox("Do Sphere-Sphere Collisions", &DoSSCollisions);
				ImGui::Checkbox("Do Sphere-Container Collisions", &DoContainerCollisions);

				PhysicsStep = ImGui::Button("Step Simulation");

				ImGui::NewLine();
				ImGui::SliderFloat("Mouse Ripple Size", &MouseRippleSize, 0.0f, 0.08f);
				ImGui::NewLine();
				ImGui::SliderInt("Substeps", &Substeps, 1, 100);
				ImGui::SliderFloat("c", &c, 0.0f, 24.0f);
				ImGui::SliderFloat("s", &s, 0.1f, 10.0f);
				ImGui::SliderFloat("Sphere Friction Coeff", &SphereFrictionCoefficient, 0.0f, 1.0f);
				ImGui::SliderFloat("Damping Coeff", &DampingCoeff, 0.01f, 4.0f);
				ImGui::SliderFloat("Alpha (Coefficient of Object-Water Interaction)", &AlphaO, 0.0f, 1.0f);
				ImGui::SliderInt("Smooth Interaction Map", &SMOOTHING_ITR, 0, 4);

				ImGui::NewLine();

				if (ImGui::Button("*Create Water Disturbance")) {
					Heightmap[int(RandomGen.Float() * Resolution * Resolution)] = 1.5;
					Heightmap[int(RandomGen.Float() * Resolution * Resolution)] = 0.5f;
				}

				if (ImGui::Button("Reset Water Sim")) {

					for (int i = 0; i < Resolution * Resolution; i++) {
						Heightmap[i] = 1.f;
					}

					memset(WaterAccelerations, 0, Resolution * Resolution * sizeof(float));
					memset(WaterVelocities, 0, Resolution * Resolution * sizeof(float));
					memset(ObjectHeights[0], 0, Resolution * Resolution * sizeof(float));
					memset(ObjectHeights[1], 0, Resolution * Resolution * sizeof(float));
				}

				ImGui::NewLine();
				ImGui::NewLine();

				ImGui::Checkbox("Wireframe", &WireFrame);
				//ImGui::Checkbox("Do Sphere Collisions", &DoCollisions);
				

				ImGui::NewLine();
				ImGui::NewLine();
				
				if (Spheres.size() < MAX_SPHERES) {
					if (ImGui::Button("PLACE Sphere")) {
						Sphere s1 = { glm::vec3(Camera.GetPosition() + Camera.GetFront() * (SphereRad * 2.1f)), glm::vec3(0.0f), CurrDens, SphereRad, glm::vec3(0.0f), {} };
						Spheres.push_back(s1);
					}

					ImGui::NewLine();

					if (ImGui::Button("THROW Sphere")) {
						Sphere s1 = { glm::vec3(Camera.GetPosition() + Camera.GetFront() * (SphereRad * 2.1f)), Camera.GetFront() * SphereVel, CurrDens, SphereRad, glm::vec3(0.0f), {}};
						Spheres.push_back(s1);
					}

					ImGui::NewLine();

					ImGui::SliderFloat("Radii of Sphere", &SphereRad, 0.05f, 2.0f);
					ImGui::SliderFloat("Sphere Density (Density of water is 1000 kg/m3)", &CurrDens, 2.0f, 4000.0f);
					ImGui::SliderFloat("Magnitude of velocity of sphere", &SphereVel, 0.0f, 16.0f);
					
				}
				ImGui::NewLine();

				if (Spheres.size() > 0) {
					if (ImGui::Button("Delete Sphere")) {
						Spheres.pop_back();
					}
				}

				ImGui::NewLine();

				if (ImGui::Button("Delete ALL Spheres")) {
					Spheres.clear();
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

			if (e.type == Simulation::EventTypes::KeyPress && e.key == GLFW_KEY_Q && this->GetCurrentFrame() > 5 && Spheres.size() < MAX_SPHERES)
			{
				Sphere s1 = { glm::vec3(Camera.GetPosition() + Camera.GetFront() * (SphereRad * 2.1f)), Camera.GetFront() * SphereVel, CurrDens, SphereRad, glm::vec3(0.0f), {} };
				Spheres.push_back(s1);
			}

			if (e.type == Simulation::EventTypes::KeyPress && e.key == GLFW_KEY_E && this->GetCurrentFrame() > 5 && Spheres.size() < MAX_SPHERES)
			{
				Sphere s1 = { glm::vec3(Camera.GetPosition() + Camera.GetFront() * (SphereRad * 2.1f)), Camera.GetFront() * 0.f, CurrDens, SphereRad, glm::vec3(0.0f), {} };
				Spheres.push_back(s1);
			}


		}


	};

	float SampleHeight(int x, int y) {
		return Heightmap[To1DIdx(x, y)];
	}

	glm::vec2 ConvertToWorldSpace(const glm::ivec2& Texel) {
		return ((glm::vec2(Texel) / float(Resolution)) * 2.0f - 1.0f) * Range;
	}

	glm::vec2 ConvertToWorldSpaceTexelCenter(const glm::ivec2& Texel) {
		return (((glm::vec2(Texel) + 0.5f) / float(Resolution)) * 2.0f - 1.0f) * Range;
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

		float Max = 0.0f;

		for (int x = 0; x < Resolution; x++) {
			for (int y = 0; y < Resolution; y++) {

				int idx = To1DIdx(x, y);

				float Velocity = WaterVelocities[To1DIdx(x, y)];
				Heightmap[idx] += Velocity * Dt;
				float DeltaH = ObjectHeights[CheckerStep][idx] - ObjectHeights[(int)(!((bool)CheckerStep))][idx];
			
				// Fix Spikes of instability 
				Max = glm::max(Max, DeltaH);

				DeltaH = glm::min(DeltaH, 0.1f);
				
				Heightmap[idx] += AlphaO * DeltaH;

				Heightmap[idx] = glm::clamp(Heightmap[idx], -float(MAX_WAVE_HEIGHT), float(MAX_WAVE_HEIGHT));
			}
		}

		//std::cout << "\n" << Max;

	}

	void SmoothObjectMap(float Dt) {


		for (int k = 0; k < SMOOTHING_ITR; k++) {
			for (int x = 0; x < Resolution; x++) {
				for (int y = 0; y < Resolution; y++) {

					float Sum = 0.0f;

					for (int p = -1; p <= 1; p++) {
						for (int q = -1; q <= 1; q++) {
							int idx = To1DIdxSafe(x + p, y + q);
							Sum += ObjectHeights[CheckerStep][idx];
						}
					}

					Sum /= 9.0f;
					ObjectHeights[CheckerStep][To1DIdxSafe(x, y)] = Sum;
				}
			}
		}

	}

	// r -> [0,1]
	float OMapWt(float r) {
		return glm::clamp(exp(-4.0f * r * r) - 0.01f, 0.0f, 1.0f);
	}

	void GenerateObjectMap(float Dt) {
		for (int x = 0; x < Resolution; x++) {
			for (int y = 0; y < Resolution; y++) {
				int idx = To1DIdx(x, y);

				// World space position of grid cell center
				glm::vec2 WorldSpace = ConvertToWorldSpaceTexelCenter(glm::ivec2(x, y));

				float NetVolumeDisplaced = 0.0f;

				for (int i = 0; i < Spheres.size(); i++) {
					glm::vec3 rsi = RSI(glm::vec3(WorldSpace.x, 1.0f, WorldSpace.y) - Spheres[i].Position, glm::vec3(0.0f, -1.0f, 0.0f), Spheres[i].Radius);
					
					float h = rsi.z;


					if (rsi.x < 0.0f && rsi.y < 0.0f) {
						h = 0.0f;
					}
					else if (rsi.x > 0.0f && rsi.y > 0.0f) {
						h = std::abs(rsi.x - rsi.y);
					}

					//h += Heightmap[idx] - 1.0f;

					h = glm::max(h, 0.0f);


					float Volume = ColumnWidth * ColumnWidth * h;
					float Fb = RhoWater * 9.81f * Volume;
					Spheres[i].NetForce += (Fb * glm::vec3(0.0f, 1.0f, 0.0f)); // Upward buoyant force
					NetVolumeDisplaced += h * OMapWt(glm::clamp(rsi.z / 8.0f, 0.0f, 1.0f));;
				}

				(ObjectHeights[CheckerStep])[idx] = NetVolumeDisplaced;
			}
		}
	}

	void CollideObjects(float dt) {


		for (int i = 0; i < Spheres.size(); i++) {

			auto& e = Spheres[i];

			if (DoSSCollisions) {
				for (int j = 0; j < Spheres.size(); j++) {

					if (i == j) {
						continue;
					}

					auto& e2 = Spheres[j];

					glm::vec3 DeltaP = e.Position - e2.Position;

					float Length = glm::length(DeltaP);

					glm::vec3 Dir = DeltaP / glm::max(Length, 0.00001f);

					if (Length <= e.Radius + e2.Radius) {

						float Delta = (e.Radius + e2.Radius) - Length;

						glm::vec3 p1, p2;

						p1 = e.Position;
						p2 = e2.Position;

						e.Position += Delta * Dir * 0.5f;
						e2.Position -= Delta * Dir * 0.5f;

					}


				}
			}

			{
				float PlayerRadius = 0.08f;

				glm::vec3 DeltaP = e.Position - Camera.GetPosition();

				float Length = glm::length(DeltaP);

				// DESPAWN IF TOO FAR
				if (Length > 4.0f * Range) {
					Spheres.erase(Spheres.begin() + i);
				}

				else {

					glm::vec3 Dir = DeltaP / glm::max(Length, 0.00001f);

					if (Length <= e.Radius + PlayerRadius) {

						float Delta = (e.Radius + PlayerRadius) - Length;

						glm::vec3 p1, p2;

						p1 = e.Position;

						e.Position += Delta * Dir;

					}
				}
			}

			if (DoContainerCollisions)
			{
				float Length = glm::length(e.Position);
				glm::vec3 Dir = e.Position / glm::max(Length, 0.00001f);

				glm::vec3 p = e.Position;

				if (Length >= 16. - e.Radius) {

					float Delta = (16. - e.Radius) - Length;
					e.Position += Delta * Dir;
					e.Velocity = (e.Position - p) / dt;
				}
			}



		}

	}


	void SimulateObjects(float dt) {

		const float g = 9.81;
		const glm::vec3 jcap = glm::vec3(0.0f, 1.0f, 0.0f);

		for (int i = 0; i < Spheres.size(); i++) {

			Sphere& s = Spheres[i];

			s.NetForce += (s.Mass() * g * -jcap);

			// Energy Loss 
			glm::vec3 OpposingDirection = - (s.Velocity / glm::max(glm::length(s.Velocity), 0.001f));
			s.NetForce += SphereFrictionCoefficient * s.Mass() * 9.81f * OpposingDirection;

			s.NetAcceleration = s.NetForce / s.Mass();
			
			s.NetForce = glm::vec3(0.0f);
			
			s.Velocity += s.NetAcceleration * dt;
			
			glm::vec3 PrevPosition = s.Position;
			
			s.Position += s.Velocity * dt;
			s.Velocity = (s.Position - PrevPosition) / dt;
		}

	}


	void DriftSpheres(float Dt) {

		for (int i = 0; i < Spheres.size(); i++) {

			Sphere& s = Spheres[i];

			for (int x = 0; x < 32; x++) {

				glm::vec2 Displaced = glm::vec2(s.Position.x, s.Position.z) + (PoissonDisk[x] * s.Radius);
				glm::vec2 UV = glm::fract(((Range + Displaced) / Range) * 0.5f);
				glm::ivec2 Texel = glm::ivec2(UV * float(Resolution));

				float Acceleration = WaterAccelerations[To1DIdxSafe(Texel.x, Texel.y)];
				float H = Heightmap[To1DIdxSafe(Texel.x, Texel.y)];

				glm::vec3 Vector = glm::vec3(s.Position - glm::vec3(Displaced.x, H, Displaced.y));

				float Length = glm::length(Vector);

				Vector /= glm::max(Length, 0.00001f);

				float Weight = Length / s.Radius;

				s.NetForce += Weight * glm::vec3(0.,1.,0.) * s.Mass() * Vector * Acceleration * DebugVar;
			}

		}

	}


	void SimulateWater(float Dt) {

		float Ddt = Dt / float(Substeps);

		for (int i = 0; i < Substeps; i++) {

			GenerateObjectMap(Ddt);
			SmoothObjectMap(Ddt);
			SimulateWaterAcceleration(Ddt);
			SimulateWaterVelocities(Ddt);
			SimulateWaterHeights(Ddt);
			//DriftSpheres(Ddt);

			if (DoContainerCollisions || DoSSCollisions)
				CollideObjects(Ddt);

			SimulateObjects(Ddt);
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

		GLClasses::CubeTextureMap Skybox;
		Skybox.CreateCubeTextureMap({
		   "Res/Sky/right.bmp",
		   "Res/Sky/left.bmp",
		   "Res/Sky/top.bmp",
		   "Res/Sky/bottom.bmp",
		   "Res/Sky/front.bmp",
		   "Res/Sky/back.bmp"
			});

		GLClasses::Texture PoolTexture;
		PoolTexture.CreateTexture("Res/pool.jpg",false, false, false, GL_TEXTURE_2D, GL_LINEAR, GL_LINEAR);

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

					float Cx = (1. / (float)(Resolution)) * float(x);
					float Cy = (1. / (float)(Resolution)) * float(y);

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
		GLClasses::Shader& RaytracingShader = ShaderManager::GetShader("RTSHADER");
		GLClasses::Shader& PostShader = ShaderManager::GetShader("POST");
		
		GLClasses::Framebuffer GBuffer[2] = { GLClasses::Framebuffer(16, 16, {{GL_RGBA16F, GL_RGBA, GL_FLOAT, false, false}}, true, true), GLClasses::Framebuffer(16, 16, {{GL_RGBA16F, GL_RGBA, GL_FLOAT, false, false}}, true, true) };

		// Spheres
		GLuint SphereSSBO = 0;
		glGenBuffers(1, &SphereSSBO);
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, SphereSSBO);
		glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(RenderSphere) * MAX_SPHERES, (void*)0, GL_DYNAMIC_DRAW);
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

		GLuint FocusSSBO = 0;
		glGenBuffers(1, &FocusSSBO);
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, FocusSSBO);
		glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(glm::vec4) * 1, (void*)0, GL_DYNAMIC_DRAW);
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

		// Create Heightmaps
		Heightmap = new float[Resolution * Resolution];
		ObjectHeights[0] = new float[Resolution * Resolution];
		ObjectHeights[1] = new float[Resolution * Resolution];
		WaterVelocities = new float[Resolution * Resolution];
		WaterAccelerations = new float[Resolution * Resolution];

		memset(Heightmap, 0, Resolution * Resolution * sizeof(float));
		memset(WaterAccelerations, 0, Resolution * Resolution * sizeof(float));
		memset(WaterVelocities, 0, Resolution * Resolution * sizeof(float));
		memset(ObjectHeights[0], 0, Resolution * Resolution * sizeof(float));
		memset(ObjectHeights[1], 0, Resolution * Resolution * sizeof(float));

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

		Sphere s1 = { glm::vec3(0.0f), glm::vec3(0.0f), 1.0f, 0.5f, glm::vec3(0.0f), {} };
		//Spheres.push_back(s1);

		glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
		glClear(GL_COLOR_BUFFER_BIT);

		SimulateObjects(0.00001f);
		SimulateWater(0.00001f);

		while (!glfwWindowShouldClose(app.GetWindow())) {

			// Update 
			glDisable(GL_CULL_FACE);

			CheckerStep = app.GetCurrentFrame() % 2;
			kProportionality = (c * c) / (s * s);
			ColumnWidth = (2.0f * Range) / float(Resolution);
			app.OnUpdate();

			// Mouse Ripple Update
			glm::vec4 FocusedUV;
			glBindBuffer(GL_SHADER_STORAGE_BUFFER, FocusSSBO);
			glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, sizeof(glm::vec4), &FocusedUV);;
			glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

			if (glfwGetMouseButton(app.GetWindow(), GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS && app.GetCurrentFrame() > 100) {
				glm::ivec2 px = glm::ivec2(glm::floor(glm::vec2(FocusedUV.x, FocusedUV.y) * glm::vec2(Resolution)));
				Heightmap[To1DIdxSafe(px.x, px.y)] += MouseRippleSize;
			}

			// FBO Update
			GBuffer[0].SetSize(app.GetWidth(), app.GetHeight());
			GBuffer[1].SetSize(app.GetWidth(), app.GetHeight());

			// Player
			MainPlayer.OnUpdate(app.GetWindow(), DeltaTime, MainPlayer.Speed, app.GetCurrentFrame());

			// SIMULATE
			if (DoSim || PhysicsStep)
			{
				SimulateWater(DeltaTime);
				PhysicsStep = false;
			}


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
			glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * Resolution * Resolution, Heightmap, GL_DYNAMIC_DRAW);

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
			RaytracingShader.Use();

			RaytracingShader.SetInteger("u_Texture", 0);
			RaytracingShader.SetInteger("u_Depth", 1);
			RaytracingShader.SetInteger("u_Skybox", 2);
			RaytracingShader.SetInteger("u_PoolTexture", 3);
			
			RaytracingShader.SetInteger("u_Spheres", Spheres.size());
			RaytracingShader.SetBool("u_RenderSpheres", RenderSpheres);
			RaytracingShader.SetBool("u_RenderPool", RenderPool);

			RaytracingShader.SetFloat("u_zNear", Camera.GetNearPlane());
			RaytracingShader.SetFloat("u_zFar", Camera.GetFarPlane());

			RaytracingShader.SetFloat("u_PoolRange", PoolRange);
			RaytracingShader.SetFloat("u_PoolHeight", PoolHeight);
			RaytracingShader.SetFloat("u_Range", Range);

			RaytracingShader.SetFloat("u_WaterBlueness", WaterBlueness);

			RaytracingShader.SetMatrix4("u_InverseProjection", glm::inverse(Camera.GetProjectionMatrix()));
			RaytracingShader.SetMatrix4("u_InverseView", glm::inverse(Camera.GetViewMatrix()));

			RaytracingShader.SetInteger("u_ScreenResH", app.GetHeight());
			RaytracingShader.SetInteger("u_ScreenResW", app.GetWidth());

			if (app.GetCursorLocked()) {
				RaytracingShader.SetInteger("u_MouseY", app.GetHeight() / 2);
				RaytracingShader.SetInteger("u_MouseX", app.GetWidth() / 2);
			}

			else {
				double mxx, myy;
				glfwGetCursorPos(app.GetWindow(), &mxx, &myy);
				myy = (double)app.GetHeight() - myy;
				glm::ivec2 fp = glm::vec2((float)mxx, (float)myy);

				RaytracingShader.SetInteger("u_MouseY", fp.y);
				RaytracingShader.SetInteger("u_MouseX", fp.x);
			}
			
			glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, SphereSSBO);
			glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, FocusSSBO);
			
			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, GBuffer[0].GetTexture());

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, GBuffer[0].GetDepthBuffer());

			glActiveTexture(GL_TEXTURE2);
			glBindTexture(GL_TEXTURE_CUBE_MAP, Skybox.GetID());

			glActiveTexture(GL_TEXTURE3);
			glBindTexture(GL_TEXTURE_2D, PoolTexture.GetTextureID());
			
			ScreenQuadVAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			ScreenQuadVAO.Unbind();

			// Blit

			glBindFramebuffer(GL_FRAMEBUFFER, 0);

			glPolygonMode(GL_FRONT, GL_FILL);
			glPolygonMode(GL_BACK, GL_FILL);

			PostShader.Use();

			PostShader.SetInteger("u_Texture", 0);

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