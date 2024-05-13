#pragma once 

#include <iostream>
#include <glm/glm.hpp>

struct Object {
	glm::vec4 Position; // w component has radius 
	glm::vec4 Velocity; 
	glm::vec4 Acceleration;
};