#version 430 core 

layout (location = 0) out vec4 o_Color;

in vec2 v_TexCoords;

uniform sampler2D u_Texture;
uniform sampler2D u_Depth;

uniform float u_zNear;
uniform float u_zFar;

uniform int u_Spheres;

uniform mat4 u_InverseProjection;
uniform mat4 u_InverseView;

layout (std430, binding = 0) buffer SSBO_Spheres {
	vec4 SphereData[];
};

vec3 SkyColour(vec3 ray)
{
    return mix(vec3(0.7), vec3(0.0f), exp(-(ray.y + 0.1)));
}

vec3 SampleIncidentRayDirection(vec2 screenspace)
{
	vec4 clip = vec4(screenspace * 2.0f - 1.0f, -1.0, 1.0);
	vec4 eye = vec4(vec2(u_InverseProjection * clip), -1.0, 0.0);
	return normalize(vec3(u_InverseView * eye));
}

float LinearizeDepth(float depth)
{
	return (2.0 * u_zNear) / (u_zFar + u_zNear - depth * (u_zFar - u_zNear));
}

float TraceSphere(vec3 Origin, vec3 Dir, float Radius)
{
    float VoV = dot(Dir, Dir);
    float Acc = VoV * Radius * Radius;

    // Solve quadratic 
    Acc += 2.0 * Origin.x * dot(Origin.yz, Dir.yz) * Dir.x;
    Acc += 2.0 * Origin.y * Origin.z * Dir.y * Dir.z;
    Acc -= dot(Origin * Origin, vec3(dot(Dir.yz, Dir.yz), dot(Dir.xz, Dir.xz), dot(Dir.xy, Dir.xy)));
    
    // No intersect 
    if (Acc < 0.0)
    {
        return -1.0;
    }
    
    Acc = sqrt(Acc);
    
    float Dist1 = (Acc - dot(Origin, Dir)) / VoV;
    float Dist2 = -(Acc + dot(Origin, Dir)) / VoV;

    if (Dist1 >= 0.0 && Dist2 >= 0.0)
    {
        return min(Dist1, Dist2);
    }
    else
    {
        return max(Dist1, Dist2);
    }
}

vec3 WorldPosFromDepth(float depth, vec2 txc)
{
    float z = depth * 2.0 - 1.0;
    vec4 ClipSpacePosition = vec4(txc * 2.0 - 1.0, z, 1.0);
    vec4 ViewSpacePosition = u_InverseProjection * ClipSpacePosition;
    ViewSpacePosition /= ViewSpacePosition.w;
    vec4 WorldPos = u_InverseView * ViewSpacePosition;
    return WorldPos.xyz;
}

void main() {

	vec3 RayOrigin = u_InverseView[3].xyz;
	vec3 RayDirection = SampleIncidentRayDirection(v_TexCoords);

    float Depth = texture(u_Depth, v_TexCoords).x;
    float LinearDepth = LinearizeDepth(Depth);

    vec3 WorldPos = WorldPosFromDepth(Depth, v_TexCoords);
    float Dist = distance(RayOrigin, WorldPos);

	o_Color = texture(u_Texture, v_TexCoords);

    bool Transparent = false;

    if (Depth > 1.0f - 0.00001f) {
        o_Color = vec4(SkyColour(RayDirection), 1.);
    }

   for (int i = 0 ; i < u_Spheres ; i ++) {
       float T = TraceSphere(RayOrigin - SphereData[i].xyz, RayDirection, SphereData[i].w);
   
       if (T > 0.0f && T < Dist && (int(gl_FragCoord.x + gl_FragCoord.y) % 2 == 0 || !Transparent)) {
           o_Color = vec4(vec3(1.,0.,0.),1.);
           Dist = T;
       }
   }

}