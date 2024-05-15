#version 330 core

layout (location = 0) in vec2 a_Position;
layout (location = 1) in vec2 a_TexCoords;

out vec2 v_TexCoords;

uniform mat4 u_ViewProj;

void main()
{
	gl_Position = u_ViewProj * 100.0f * vec4(a_Position, 1.0f, 1.0f);
	v_TexCoords = a_TexCoords;
}