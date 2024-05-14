#pragma once 

#include <iostream>
#include <glm/glm.hpp>

struct Object {
	glm::vec2 Position;
	glm::vec2 Velocity; 
	glm::vec2 Force;
	glm::vec2 MassRadius;
	glm::vec2 Dx;
	glm::vec2 Dv;
	
};