#version 330 core

layout (location = 0) in vec2 a_Position;

uniform mat4 u_ViewProj;
uniform mat4 u_ViewProjRot;

void main()
{
	gl_Position = u_ViewProjRot * vec4(a_Position, 1.0f, 1.0f);
}