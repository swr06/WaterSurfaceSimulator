#version 330 core 

layout (location = 0) out vec4 o_Color;

uniform sampler2D u_Texture;

void main() {

	o_Color = vec4(1.,1.,0.,1.);// texture(u_Texture, v_TexCoords);

}