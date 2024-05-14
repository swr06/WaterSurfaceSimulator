#pragma once

#include <glm/glm.hpp>
#include "FpsCamera.h"

#include <GLFW/glfw3.h>

namespace Simulation
{
	class Player
	{
	public :

		Player();
		void OnUpdate(GLFWwindow* window, float dt, float speed, int frame);

		void TestCollision(glm::vec3& position, glm::vec3 vel);
		void Jump();

		FPSCamera Camera;
		bool Freefly = false;
		float Sensitivity = 0.25;
		float Speed = 0.1250f;

		glm::vec3 m_Position;
		glm::vec3 m_Velocity;
		glm::vec3 m_Acceleration;
		bool m_isOnGround;
		bool DisableCollisions = false;


	private :

	};
}