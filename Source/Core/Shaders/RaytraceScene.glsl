#version 450 core 

layout (location = 0) out vec4 o_Color;

in vec2 v_TexCoords;

uniform sampler2D u_Texture;
uniform sampler2D u_Depth;
uniform sampler2D u_PoolTexture;
uniform samplerCube u_Skybox;

uniform float u_zNear;
uniform float u_zFar;

uniform float u_PoolRange;
uniform float u_PoolHeight;

uniform float u_Range;
uniform int u_Res;

uniform int u_Spheres;

uniform bool u_RenderSpheres;
uniform bool u_RenderPool;

uniform mat4 u_InverseProjection;
uniform mat4 u_InverseView;

uniform int u_ScreenResW;
uniform int u_ScreenResH;
uniform int u_MouseX;
uniform int u_MouseY;

uniform vec3 u_SunDirection;

uniform float u_WaterBlueness;

layout (std430, binding = 0) buffer SSBO_Spheres {
	vec4 SphereData[];
};

layout (std430, binding = 1) buffer SSBO_CenterPx {
	vec4 CenterPxPacked;
};


layout (std430, binding = 2) buffer SSBO_HM {
	float Heightmap[];
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

vec3 SampleSkyAt(vec3 D) {
    return pow(texture(u_Skybox, D).xyz, vec3(1.4f));
}


int To1DIdx(int x, int y) {
	return (y * u_Res) + x;
}

float Sample(vec2 UV) {
	UV = clamp(UV, 0., 1.);
	ivec2 Texel = ivec2(UV * vec2(float(u_Res)));
	return Heightmap[To1DIdx(Texel.x, Texel.y)];
}

float SamplePx(ivec2 px) {
	return Heightmap[To1DIdx(px.x, px.y)];
}

vec2 SpiralPoint(float angle, float scale) {
    return vec2(sin(angle), cos(angle)) * pow(angle / scale, 1.0 / (sqrt(5.0) * 0.5 + 0.5));
}

float Bilinear(vec2 SampleUV)
{
    const ivec2 Cross[4] = ivec2[4](ivec2(0, 0), ivec2(1, 0), ivec2(0, 1), ivec2(1, 1));

    // Relative to center
    vec2 SamplingFragment = (SampleUV * vec2(u_Res)) - vec2(0.5f); 

    // Find how much you need to interpolate across both axis
    vec2 f = fract(SamplingFragment);

    // Fetch 4 neighbours
    float Fetch[4];
    Fetch[0] = SamplePx(ivec2(SamplingFragment) + Cross[0]);
    Fetch[1] = SamplePx(ivec2(SamplingFragment) + Cross[1]);
    Fetch[2] = SamplePx(ivec2(SamplingFragment) + Cross[2]);
    Fetch[3] = SamplePx(ivec2(SamplingFragment) + Cross[3]);

    // Interpolate first based on x position
    float Interp1 = mix(Fetch[0], Fetch[1], float(f.x));
    float Interp2 = mix(Fetch[2], Fetch[3], float(f.x));

    return mix(Interp1, Interp2, float(f.y));
}

float GetHeightAt(vec3 W) {
    vec2 UV = fract(((u_Range + W.xz) / u_Range) * 0.5f);
    return Sample(UV);
}

float GetHeightAtSmooth(vec3 W) {
    vec2 UV = fract(((u_Range + W.xz) / u_Range) * 0.5f);
    return Bilinear(UV);
}

vec3 SampleNormals(vec3 WorldPosition) {

    float Width = (u_Range / float(u_Res));

    vec3 WestVertex = WorldPosition + vec3(Width, 0.0f, 0.0f);
    WestVertex.y = GetHeightAtSmooth(WestVertex);

    vec3 NorthVertex = WorldPosition + vec3(0.0f, 0.0f, Width);
    NorthVertex.y = GetHeightAtSmooth(NorthVertex);
    return normalize(cross(WorldPosition - NorthVertex, WorldPosition - WestVertex));
}

float CalculateCausticsNV(in vec3 worldPosition, in float depth) 
{
    const float PI = 3.14159f;
    float samples = 16.0;
    float focus = 0.9f;
    float surfaceDistanceUp = depth * abs(u_SunDirection.y);

    float filterRadius = 0.2;

    float radius = filterRadius * depth;
    float inverseDistanceThreshold = sqrt(samples / PI) * focus / radius;

    vec3 flatRefractVector = refract(u_SunDirection, vec3(0, 1, 0), 1.00028/1.333);
    vec3 flatRefract = flatRefractVector * surfaceDistanceUp / abs(flatRefractVector.y);
    vec3 surfacePosition = worldPosition - flatRefract;

    float finalCaustic = 0.0;

    for(float i = 0; i <= samples; i++) 
    {
        vec3 samplePos = surfacePosition;
        samplePos.xz += SpiralPoint(i + 0.5f, samples) * radius;
        vec3 normal = normalize(SampleNormals(vec3(samplePos.x,1.,samplePos.z)).xyz);
        vec3 refractVector = refract(u_SunDirection, normal, 1.00028/1.333);
        samplePos = refractVector * (surfaceDistanceUp / abs(refractVector.y)) + samplePos;
        finalCaustic += clamp(1.0 - distance(worldPosition, samplePos) * inverseDistanceThreshold,0.,1.);
    }

    finalCaustic *= focus * focus;
    return pow(finalCaustic, 1.0);
}


// Based on the blog post by inigo quilez 
vec2 IntersectBox(in vec3 ro, in vec3 invrd, in vec3 rad, out vec3 oN) 
{
    vec3 m = invrd;
    vec3 n = m * ro;
    vec3 k = abs(m) * rad;
    vec3 t1 = -n - k;
    vec3 t2 = -n + k;

    float tN = max(max(t1.x, t1.y), t1.z);
    float tF = min(min(t2.x, t2.y), t2.z);
	
    if(tN > tF || tF < 0.0) return vec2(-1.0); // no intersection
    
    oN = -sign(invrd) * step(t1.yzx,t1.xyz) * step(t1.zxy,t1.xyz);

    return vec2( tN, tF );
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

float TraceSpheres(vec3 O, vec3 D) {

    if (!u_RenderSpheres) {
        return -1.;
    }

    float Dist = 100000.;
    bool Int = false;

    for (int i = 0 ; i < u_Spheres ; i ++) {
           float T = TraceSphere(O - SphereData[i].xyz, D, SphereData[i].w);
   
           if (T > 0.0f && T < Dist) {
               Dist = T;
               Int = true;
           }
       }

    return Int ? Dist : -1.;

}

// Filter with AA 
vec4 SmoothFilter(sampler2D tex, vec2 uv)
{
    vec2 res = textureSize(tex,0).xy;
    const float width = 1.;
    uv = uv * res;
    vec2 uv_floor = floor(uv + 0.5);
    vec2 uv_fract = fract(uv + 0.5);
    vec2 uv_aa = fwidth(uv) * width * 0.5;
    uv_fract = smoothstep(vec2(0.5) - uv_aa, vec2(0.5) + uv_aa,uv_fract);
    vec2 NewUV = (uv_floor + uv_fract - 0.5) / res;
    return texture(tex,NewUV);
}

// Wanted a clean way to get box uvs, i guess this works :) 
vec2 GetBoxUV(vec3 p, in vec3 n)
{
    vec3 Mapping = abs(n);
    Mapping = (max(Mapping, 0.00001));
    float sigma = (Mapping.x + Mapping.y + Mapping.z);
    Mapping /= sigma;
    return mod(p.yz * Mapping.x + p.xz * Mapping.y + p.xy * Mapping.z, 1.);
}

vec4 IntersectPool(in vec3 o, in vec3 inv_d) {
    
    if (!u_RenderPool) {
        return vec4(-1.);
    }

    const float Bias = 0.0f;

    float Thickness = 0.04f;

    const float Size = u_PoolRange;
    const float Height = u_PoolHeight;

    float c[3] = float[](1.,1.,1.);


    const vec3 BoxPositions[5] = vec3[]( vec3(0.,-0.,Size + Bias), vec3(0.,-0.,-Size - Bias), vec3(Size + Bias,-0.,0.),
                                   vec3(-Size - Bias,-0.,0.), vec3(0.,-Height,0.) ); 
    const vec3 BoxRanges[5] = vec3[]( vec3(Size, Height, Thickness), vec3(Size, Height, Thickness), vec3(Thickness, Height, Size), vec3(Thickness, Height, Size), vec3(Size, Thickness, Size) );

    float MaxDist = 100000.0f;
    vec2 intersect = vec2(-1.);
    vec3 Nf = vec3(-1.);
    bool Int = false;
    vec3 n;

    for (int i = 0 ; i < 5 ; i ++) {
        intersect = IntersectBox(o - BoxPositions[i] + vec3(0.,Height - 1.5f,0.), inv_d, BoxRanges[i], n) ;

        if (intersect.x > 0.0f && intersect.x < MaxDist) {
            MaxDist = intersect.x;
            Nf = n;
            Int = true;
        }
    }

    return Int ? vec4(MaxDist, Nf) : vec4(-1.);
}

vec3 GetPoolShading(in vec3 wp, in vec3 n) {
    vec2 Uv = GetBoxUV(wp, n);
    //float Caustic = n == vec3(0.,1.,0.) ? CalculateCausticsNV(wp, 4.0f) : 1.;
    float Caustic = 1.0f;
    vec3 Albedo = SmoothFilter(u_PoolTexture, Uv.xy).xyz;
    return Albedo * Caustic;
}

vec4 GetRayShading(in vec3 o, in vec3 dir, in vec3 invdir) {
    
    
    float SphereT = TraceSpheres(o, dir);
    vec4 Pool = IntersectPool(o, invdir);

    // both intersected 
    if (Pool.x > 0.0f && SphereT > 0.0f) {
        //if (min(SphereT, Pool.x) < PlaneT) 
        {
            bool ShadePool = Pool.x < SphereT;

            if (ShadePool) {
                return vec4(GetPoolShading(o + dir * Pool.x, Pool.yzw), Pool.x);
            }

            else {
                return vec4(vec3(1.,0.,0.), SphereT);
            }
        }
    }
    
    // Only pool
    if (Pool.x > 0.0f && SphereT < 0.0f) {
        return vec4(GetPoolShading(o + dir * Pool.x, Pool.yzw), Pool.x);
    }

    // Only sphere 
    if (Pool.x < 0.0f && SphereT > 0.0f) {
        return vec4(vec3(1.,0.,0.), SphereT);
    }

    // None
    return vec4(SampleSkyAt(dir), -1.);
}

vec3 GetPrimaryRayWater(in vec3 wp, in vec3 dir, in vec3 n) {
    vec3 ReflectedDir = reflect(dir, n);
    vec3 RefractedDir = refract(dir, n, 1.0f / 1.33f);

    vec3 Refracted = GetRayShading(wp,RefractedDir,1./RefractedDir).xyz;
    vec3 Reflected = GetRayShading(wp,ReflectedDir,1./ReflectedDir).xyz;
    float Fresnel = (pow(1.0 - max(0.0, dot(-vec3(0.,1.,0.), dir)), 2.f));;

    vec3 WaterColor = vec3(0.75f, 0.9f, 1.6f);
    //vec3 WaterColor = sqrt(texture(u_Skybox, vec3(0.,1.,0.)).xyz);
    return mix(Refracted * clamp(WaterColor / vec3(u_WaterBlueness, u_WaterBlueness,1.), 0. ,1.), Reflected, Fresnel);;
}

void GetPrimaryRayShading(in vec3 o, in vec3 dir, in vec3 invdir, inout float primtraversal, in vec3 Normals, in vec3 WP, inout vec3 oColor) {
    
    float SphereT = TraceSpheres(o, dir);
    vec4 Pool = IntersectPool(o, invdir);
    bool ShadeWater = true;

    float WorldTraversal = primtraversal < 0. ? 100000.0f : primtraversal;

    // both intersected 
    if (Pool.x > 0.0f && SphereT > 0.0f) {
        if (min(SphereT, Pool.x) < WorldTraversal) 
        {
            bool ShadePool = Pool.x < SphereT;

            if (ShadePool) {
                oColor.xyz = GetPoolShading(o + dir * Pool.x, Pool.yzw);
                primtraversal = Pool.x;
            }

            else {
                oColor = vec3(1.,0.,0.);
                primtraversal = SphereT;
            }

            ShadeWater = false;
        }
    }
    
    // Only pool
    else if (Pool.x > 0.0f && SphereT < 0.0f && Pool.x < WorldTraversal) {
        oColor.xyz = GetPoolShading(o + dir * Pool.x, Pool.yzw);
        primtraversal = Pool.x;
        ShadeWater = false;
    }

    // Only sphere 
    else if (Pool.x < 0.0f && SphereT > 0.0f && SphereT < WorldTraversal) {
        oColor = vec3(1.,0.,0.);
        primtraversal = SphereT;
        ShadeWater = false;
    }

    if (ShadeWater && primtraversal > 0.0f) {
        oColor = GetPrimaryRayWater(WP,dir,Normals);
    }

    else if (primtraversal < 0.0f) {
        oColor = pow(SampleSkyAt(dir), vec3(1.5f));
    }
}

void main() {

	vec3 RayOrigin = u_InverseView[3].xyz;
	vec3 RayDirection = (SampleIncidentRayDirection(v_TexCoords));
    vec3 InvDir = 1.0f / RayDirection;

    vec3 n = vec3(0.0f);

    float eps = 0.1f;
    float range = 4.0f;


    //

    float Depth = texture(u_Depth, v_TexCoords).x;
    float LinearDepth = LinearizeDepth(Depth);

    vec3 WorldPos = WorldPosFromDepth(Depth, v_TexCoords);

    vec3 R = normalize(WorldPos - RayOrigin);

    bool NonWaterPx = Depth > 1.0f - 0.00001f;

    // Ripple handling
     if (ivec2(gl_FragCoord.xy) == ivec2(u_MouseX, u_MouseY)) {
        vec2 wsUV = fract(((u_Range + WorldPos.xz) / u_Range) * 0.5f);
        CenterPxPacked = vec4(wsUV, float(NonWaterPx), 1.);
    }

    float Dist = NonWaterPx ? -1. : distance(RayOrigin, WorldPos);

	vec3 NormalsSampled = normalize(texture(u_Texture, v_TexCoords).xyz);

    NormalsSampled *= sign(dot(NormalsSampled, vec3(0.0f, 1.0f, 0.0f)));
    
    vec3 Color;
    GetPrimaryRayShading(RayOrigin, RayDirection, InvDir, Dist, NormalsSampled, WorldPos, Color);

    o_Color = vec4(Color, 1.);

    //if (!NonWaterPx)
    //     o_Color = vec4(CalculateCausticsNV(WorldPos.xyz, 4.0f).xxx, 1.);
}