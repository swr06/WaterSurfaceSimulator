#version 330 core 

layout (location = 0) out vec4 o_Color;

in vec2 v_TexCoords;

uniform sampler2D u_Texture;

void main() {

	o_Color = vec4(1.,0.,0.,1.);// texture(u_Texture, v_TexCoords);

}