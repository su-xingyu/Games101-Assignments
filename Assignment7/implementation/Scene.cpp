//
// Created by Göksu Güvendiren on 2019-05-14.
//

#include "Scene.hpp"


void Scene::buildBVH() {
    printf(" - Generating BVH...\n\n");
    this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
}

Intersection Scene::intersect(const Ray &ray) const
{
    return this->bvh->Intersect(ray);
}

void Scene::sampleLight(Intersection &pos, float &pdf) const
{
    float emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
        }
    }
    float p = get_random_float() * emit_area_sum;
    emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
            if (p <= emit_area_sum){
                objects[k]->Sample(pos, pdf);
                break;
            }
        }
    }
}

bool Scene::trace(
        const Ray &ray,
        const std::vector<Object*> &objects,
        float &tNear, uint32_t &index, Object **hitObject)
{
    *hitObject = nullptr;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        float tNearK = kInfinity;
        uint32_t indexK;
        Vector2f uvK;
        if (objects[k]->intersect(ray, tNearK, indexK) && tNearK < tNear) {
            *hitObject = objects[k];
            tNear = tNearK;
            index = indexK;
        }
    }


    return (*hitObject != nullptr);
}

// Implementation of Path Tracing
Vector3f Scene::castRay(const Ray &ray, int depth) const
{
    // TO DO Implement Path Tracing Algorithm here
    Intersection inter_obj = this->intersect(ray);
    if (!inter_obj.happened) {
        return Vector3f(0);
    }

    if (inter_obj.m->hasEmission()) {
        return inter_obj.m->getEmission();  // return directly since only light source in this scene emit light,
                            // otherwise should regard it as L_e
    }

    Vector3f wo = -ray.direction;
    Vector3f p = inter_obj.coords;
    Vector3f N = inter_obj.normal;

    Vector3f L_dir(0);

    Intersection inter_light;
    float pdf_light = 0.f;
    this->sampleLight(inter_light, pdf_light);

    Vector3f x = inter_light.coords;
    Vector3f ws = (x-p).normalized();
    Vector3f NN = inter_light.normal;

    Intersection inter_block = this->intersect(Ray(p, ws));

    if ((inter_block.coords-x).norm() < 0.01) {
        float cos_theta_p = std::max(0.f, dotProduct(ws, N));
        float cos_theta_x = std::max(0.f, dotProduct(-ws, NN));
        L_dir += inter_light.emit * inter_obj.m->eval(ws, wo, N) * cos_theta_p * cos_theta_x / std::pow((x-p).norm(), 2) / pdf_light;
    }

    Vector3f L_indir(0);

    if (get_random_float() < this->RussianRoulette) {
        Vector3f wi = inter_obj.m->sample(wo, N);
        L_indir += this->castRay(Ray(p, wi), depth) * inter_obj.m->eval(wi, wo, N) * std::max(0.f, dotProduct(wi, N)) / inter_obj.m->pdf(wi, wo, N) / this->RussianRoulette;
    }
    return L_dir + L_indir;
}