#include "raytracer.hpp"

#include "args.hpp"
#include "Camera.h"
#include "hit.hpp"
#include "lights.hpp"
#include "material.hpp"
#include "objects.hpp"
#include "ray.hpp"
#include "SceneParser.h"

#define EPSILON 0.001f

using namespace std;
using namespace FW;

namespace {

// Compute the mirror direction for the incoming direction, given the surface normal.
Vec3f mirrorDirection(const Vec3f& normal, const Vec3f& incoming) {
	// YOUR CODE HERE (R8)
	// Pay attention to the direction which things point towards, and that you only
	// pass in normalized vectors.

	return Vec3f(incoming - 2.0f * (incoming.dot(normal)) * normal).normalized();
}

bool transmittedDirection(const Vec3f& normal, const Vec3f& incoming, 
		float index_i, float index_t, Vec3f& transmitted) {
	// YOUR CODE HERE (EXTRA)
	// Compute the transmitted direction for the incoming direction, given the surface normal
	// and indices of refraction. Pay attention to the direction which things point towards!
	// You should only pass in normalized vectors!
	// The function should return true if the computation was successful, and false
	// if the transmitted direction can't be computed due to total internal reflection.
	return true;
}

} // namespace

Vec3f RayTracer::traceRay(Ray& ray, float tmin, int bounces, float refr_index, Hit& hit, FW::Vec3f debug_color) const {
	// initialize a hit to infinitely far away
	hit = Hit(FLT_MAX);

	// Ask the root node (the single "Group" in the scene) for an intersection.
	bool intersect = false;
	if(scene_.getGroup())
		intersect = scene_.getGroup()->intersect(ray, hit, tmin);

	// Write out the ray segment if visualizing (for debugging purposes)
	if (debug_trace)
		debug_rays.push_back(RaySegment(ray.origin, ray.direction.normalized()*min(100.0f,hit.t), hit.normal, debug_color));

	// if the ray missed, we return the background color.
	if (!intersect)
		return scene_.getBackgroundColor();
	
	Material* m = hit.material;
	assert(m != nullptr);

	// get the intersection point and normal.
	Vec3f normal = hit.normal;
	Vec3f point = ray.pointAtParameter(hit.t);

	// YOUR CODE HERE (R1)
	// Apply ambient lighting using the ambient light of the scene
	// and the diffuse color of the material.
	//Vec3f answer = Vec3f(1.0f);
	//Vec3f answer = Vec3f(intersect,scene_.getAmbientLight,m->diffuse_color);
	//Vec3f answer = m->diffuse_color(point);
	Vec3f answer = (scene_.getAmbientLight()*(m->diffuse_color(point)));

	// YOUR CODE HERE (R4 & R7)
	// For R4, loop over all the lights in the scene and add their contributions to the answer.
	// If you wish to include the shadow rays in the debug visualization, insert the segments
	// into the debug_rays list similar to what is done to the primary rays after intersection.
	// For R7, if args_.shadows is on, also shoot a shadow ray from the hit point to the light
	// to confirm it isn't blocked; if it is, ignore the contribution of the light.

	for (auto i = 0; i < scene_.getNumLights(); i++) {
		float distance=0.0f;
		Vec3f dir_to_light{ 0.0f,0.0f,0.0f }, incident_intensity{ 0.0f,0.0f,0.0f };
		scene_.getLight(i)->getIncidentIllumination(point, dir_to_light, incident_intensity, distance);
		
		if (args_.shadows) {
			Ray ray_shadow = Ray(point, dir_to_light);
			Hit hit2 = Hit(distance);
			bool shadow = scene_.getGroup()->intersect(ray_shadow, hit2, EPSILON);

			if (!shadow)
				answer += m->shade(ray, hit, dir_to_light, incident_intensity, args_.shade_back);
		}
		else
			answer += m->shade(ray, hit, dir_to_light, incident_intensity, args_.shade_back);
	}

	// are there bounces left?
	if (bounces >= 1) {
		// reflection, but only if reflective coefficient > 0!
		if (m->reflective_color(point).length() > 0.0f) {
			// YOUR CODE HERE (R8)
			// Generate and trace a reflected ray to the ideal mirror direction and add
			// the contribution to the result. Remember to modulate the returned light
			// by the reflective color of the material of the hit point.

			Ray ray_mirror = Ray(point, mirrorDirection(hit.normal, ray.direction));
			Hit hit_mirror = Hit();

			answer += m->reflective_color(point) * traceRay(ray_mirror, EPSILON, --bounces, refr_index, hit_mirror, debug_color);	
		}

		// refraction, but only if surface is transparent!
		if (m->transparent_color(point).length() > 0.0f) {
			// YOUR CODE HERE (EXTRA)
			// Generate a refracted direction and trace the ray. For this, you need
			// the index of refraction of the object. You should consider a ray going through
			// the object "against the normal" to be entering the material, and a ray going
			// through the other direction as exiting the material to vacuum (refractive index=1).
			// (Assume rays always start in vacuum, and don't worry about multiple intersecting
			// refractive objects!) Remembering this will help you figure out which way you
			// should use the material's refractive index. Remember to modulate the result
			// with the material's refractiveColor().
			// REMEMBER you need to account for the possibility of total internal reflection as well.
		}
	}
	return answer;
}
