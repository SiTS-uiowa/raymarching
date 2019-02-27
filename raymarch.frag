const int MAX_MARCHING_STEPS = 250;
const float MIN_DIST = 0.0;
const float MAX_DIST = 100.0;
const float EPSILON = 0.0001;
const float PI = 3.1415926535898;

float intersectSDF(float distA, float distB) {
    return max(distA, distB);
}

// polynomial smooth min (k = 0.1);
float smin( float a, float b, float k )
{
    float h = max( k-abs(a-b), 0.0 )/k;
    return min( a, b ) - h*h*k*(1.0/4.0);
}

// cubic smooth min (k = 0.1);
float sminCubic( float a, float b, float k )
{
    float h = max( k-abs(a-b), 0.0 )/k;
    return min( a, b ) - h*h*h*k*(1.0/6.0);
}

float unionSDF(float distA, float distB) {
    return sminCubic(distA, distB, 0.1);
}

float differenceSDF(float distA, float distB) {
    return max(distA, -distB);
}

float sphereSDF(vec3 p, float s) {
	return length(p) - s;
}

float boxSDF(vec3 p, vec3 b) {
    vec3 d = abs(p) - b;
    return length(max(d,0.0)) + min(max(d.x,max(d.y,d.z)),0.0);
}

float displacement(vec3 p) {
    return sin(20.*(p.x + iTime*0.01))*sin(20.*(p.y + iTime*0.02))*sin(20.*(p.z + iTime*0.03))*0.01;
}

float sceneSDF(vec3 samplePoint) {
    vec3 c = vec3(5.,0.,5.);
    vec3 q = mod(samplePoint,c)-0.5*c;
    
    float sphereDist = sphereSDF(q, 1.2 + sin(iTime)*0.2);
    float cubeDist = boxSDF(q, vec3(1.,1.,1.));
    float orbiterSpeed = 0.2;
    float orbiterSize = 0.25;
    float orbiter1 = sphereSDF(q + vec3(sin(iTime*orbiterSpeed)*2., 0., cos(iTime*orbiterSpeed)*2.), orbiterSize);
    float orbiter2 = sphereSDF(q + vec3(sin(iTime*orbiterSpeed+0.5*PI)*2., 0., cos(iTime*orbiterSpeed+0.5*PI)*2.), orbiterSize);
    float orbiter3 = sphereSDF(q + vec3(sin(iTime*orbiterSpeed+PI)*2., 0., cos(iTime*orbiterSpeed+PI)*2.), orbiterSize);
    float orbiter4 = sphereSDF(q + vec3(sin(iTime*orbiterSpeed+1.5*PI)*2., 0., cos(iTime*orbiterSpeed+1.5*PI)*2.), orbiterSize);

    float orbiter5 = sphereSDF(q + vec3(sin(-1.*iTime*orbiterSpeed)*2., 0., cos(-1.*iTime*orbiterSpeed)*2.), orbiterSize);
    float orbiter6 = sphereSDF(q + vec3(sin(-1.*iTime*orbiterSpeed+0.5*PI)*2., 0., cos(-1.*iTime*orbiterSpeed+0.5*PI)*2.), orbiterSize);
    float orbiter7 = sphereSDF(q + vec3(sin(-1.*iTime*orbiterSpeed+PI)*2., 0., cos(-1.*iTime*orbiterSpeed+PI)*2.), orbiterSize);
    float orbiter8 = sphereSDF(q + vec3(sin(-1.*iTime*orbiterSpeed+1.5*PI)*2., 0., cos(-1.*iTime*orbiterSpeed+1.5*PI)*2.), orbiterSize);

    float orbiters = unionSDF(orbiter1,unionSDF(orbiter2,unionSDF(orbiter3,unionSDF(orbiter4,unionSDF(orbiter5,unionSDF(orbiter6,unionSDF(orbiter7,orbiter8)))))));
    return unionSDF(unionSDF(cubeDist, sphereDist), orbiters) + displacement(samplePoint);
}

/**
 * Return the shortest distance from the eyepoint to the scene surface along
 * the marching direction. If no part of the surface is found between start and end,
 * return end.
 * 
 * eye: the eye point, acting as the origin of the ray
 * marchingDirection: the normalized direction to march in
 * start: the starting distance away from the eye
 * end: the max distance away from the ey to march before giving up
 */
float shortestDistanceToSurface(vec3 eye, vec3 marchingDirection, float start, float end) {
    float depth = start;
    for (int i = 0; i < MAX_MARCHING_STEPS; i++) {
        float dist = sceneSDF(eye + depth * marchingDirection);
        if (dist < EPSILON) {
			return depth;
        }
        depth += dist;
        if (depth >= end) {
            return end;
        }
    }
    return end;
}
            

/**
 * Return the normalized direction to march in from the eye point for a single pixel.
 * 
 * fieldOfView: vertical field of view in degrees
 * size: resolution of the output image
 * fragCoord: the x,y coordinate of the pixel in the output image
 */
vec3 rayDirection(float fieldOfView, vec2 size, vec2 fragCoord) {
    vec2 xy = fragCoord - size / 2.0;
    float z = size.y / tan(radians(fieldOfView) / 2.0);
    return normalize(vec3(xy, -z));
}

vec3 estimateNormal(vec3 p) {
    return normalize(vec3(
        sceneSDF(vec3(p.x + EPSILON, p.y, p.z)) - sceneSDF(vec3(p.x - EPSILON, p.y, p.z)),
        sceneSDF(vec3(p.x, p.y + EPSILON, p.z)) - sceneSDF(vec3(p.x, p.y - EPSILON, p.z)),
        sceneSDF(vec3(p.x, p.y, p.z  + EPSILON)) - sceneSDF(vec3(p.x, p.y, p.z - EPSILON))
    ));
}

mat3 viewMatrix(vec3 eye, vec3 center, vec3 up) {
    // Based on gluLookAt man page
    vec3 f = normalize(center - eye);
    vec3 s = normalize(cross(f, up));
    vec3 u = cross(s, f);
    return mat3(s, u, -f);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 aspect = vec2(iResolution.x/iResolution.y, 1.0);
	vec2 screenCoords = (2.0*fragCoord.xy/iResolution.xy - 1.0)*aspect;
    
    float theta = iTime * 0.1;
    float camDist = 4.;
    
    vec3 rd = rayDirection(90.0, iResolution.xy, fragCoord);
    vec3 ro = vec3(camDist * sin(theta), 2., camDist * cos(theta));
    
    mat3 viewToWorld = viewMatrix(ro, vec3(0., 0., 0.), vec3(0.0, 1.0, 0.0));
    
    rd = viewToWorld * rd;
    
    vec3 bgcolor = vec3(1.,0.97,0.92)*0.15;
    float bgshade = (1.0-length(vec2(screenCoords.x/aspect.x, screenCoords.y+0.5) )*0.8);
	bgcolor *= bgshade;
    
    float dist = shortestDistanceToSurface(ro, rd, MIN_DIST, MAX_DIST);
    if ( dist >= MAX_DIST ) {
	    fragColor = vec4(bgcolor, 1.0);
	    return;
	}
    /*
    vec3 sp = ro + rd*dist;
    vec3 surfNormal = estimateNormal(sp);
    fragColor = vec4(surfNormal * 0.5 + vec3(0.5,0.5,0.5), 1.0);
	*/
    
    vec3 sp = ro + rd*dist;
    vec3 surfNormal = estimateNormal(sp);
    
    vec3 lp = vec3(sin(theta) * 1.5, 2., cos(theta) * 1.5);
    vec3 ld = lp-sp;
    vec3 lcolor = vec3(1.,0.97,0.92);
    
    float len = length( ld ); // Distance from the light to the surface point.
	ld /= len; // Normalizing the light-to-surface, aka light-direction, vector.
	float lightAtten = min( 1.0 / ( 0.25*len*len ), 1.0 ); // Keeps things between 0 and 1.
    
    vec3 ref = reflect(-ld, surfNormal); 
    vec3 sceneColor = vec3(0.0);
    vec3 objColor = vec3(0.7, 0.3, 1.0);
    
    float ambient = .1; //The object's ambient property. You can also have a global and light ambient property, but we'll try to keep things simple.
	float specularPower = 16.0; // The power of the specularity. Higher numbers can give the object a harder, shinier look.
	float diffuse = max( 0.0, dot(surfNormal, ld) ); //The object's diffuse value, which depends on the angle that the light hits the object.
	//The object's specular value, which depends on the angle that the reflected light hits the object, and the viewing angle... kind of.
	float specular = max( 0.0, dot( ref, normalize(ro-sp)) ); 
	specular = pow(specular, specularPower); // Ramping up the specular value to the specular power for a bit of shininess.
		
	// Bringing all the lighting components together to color the screen pixel. By the way, this is a very simplified version of Phong lighting. 
	// It's "kind of" correct, and will suffice for this example. After all, a lot of lighting is fake anyway.
	sceneColor += (objColor*(diffuse*0.8+ambient)+specular*0.5)*lcolor*lightAtten;
    
    fragColor = vec4(clamp(sceneColor, 0.0, 1.0), 1.0);
}
