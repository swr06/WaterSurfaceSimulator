#version 450 core 

#define PI 3.14159265359 

layout(local_size_x = 16, local_size_y = 1, local_size_z = 1) in;
layout(rgba16f, binding = 0) uniform image2D o_OutputData;

uniform float u_Dt;
uniform float u_Time;
uniform int u_ObjectCount;

struct Object {
	vec4 Position; // w component has radius 
	vec4 Velocity; 
	vec4 Acceleration;
};

layout (std430, binding = 0) buffer ObjectSSBO {
	Object SimulationObjects[];
};

float remap(float x, float a, float b, float c, float d)
{
    return (((x - a) / (b - a)) * (d - c)) + c;
}

float HASH2SEED = 0.0f;
vec2 hash2() 
{
	return fract(sin(vec2(HASH2SEED += 0.1, HASH2SEED += 0.1)) * vec2(43758.5453123, 22578.1459123));
}

void main() {

	int Index = int(gl_GlobalInvocationID.x);

	if (Index > u_ObjectCount) {
		return;
	}

	float dt = u_Dt;

	Object CurrentObject = SimulationObjects[Index];
	Object Unupdated = SimulationObjects[Index];

	// s = s0 + ut + 1/2at^2
	CurrentObject.Position.xyz = CurrentObject.Position.xyz + CurrentObject.Velocity.xyz * dt + CurrentObject.Acceleration.xyz * dt * dt;
	CurrentObject.Acceleration = vec4(0.0f, -9.8f, 0.0f, 0.0f) * 400.0f;
	CurrentObject.Velocity.xyz = (CurrentObject.Position.xyz - Unupdated.Position.xyz) / dt;

	vec3 ToObject = CurrentObject.Position.xyz - vec3(0.0f);
	float Length = length(ToObject);
	vec3 Normal = ToObject / Length;

	if (Length > 350 - CurrentObject.Position.w) {
		//CurrentObject.Velocity.xyz = -Normal * (Length - CurrentObject.Position.w);
		CurrentObject.Velocity.xyz += -Normal * (Length - CurrentObject.Position.w);
	}

	SimulationObjects[Index] = CurrentObject;
}
