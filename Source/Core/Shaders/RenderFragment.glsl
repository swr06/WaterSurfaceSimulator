#version 330 core 

layout (location = 0) out vec4 o_Color;

in vec2 v_TexCoords;
in vec3 v_Position;
in vec3 v_RawPosition;
in float v_Radius;

uniform sampler2D u_Texture;
uniform vec2 u_Dims;

void main() {

	vec2 UV = gl_FragCoord.xy / u_Dims;

	float Aspect = u_Dims.x / u_Dims.y;
	vec2 AspectCorrect = vec2(Aspect, 1.0f);

	float Distance = distance(v_RawPosition.xy * AspectCorrect, v_Position.xy * AspectCorrect);
    float Gradient = 0.000001f; //fwidth(d);
    float Circle = smoothstep(v_Radius + Gradient, v_Radius - Gradient, Distance);

	if (Circle < 0.000001f) {
		discard;
	}

	o_Color = vec4(Circle.xxx,1.);

}