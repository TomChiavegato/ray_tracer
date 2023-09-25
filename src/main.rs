use std::fs::File;
use std::io;
use std::io::prelude::*;
use itertools::Itertools;
use indicatif::ParallelProgressIterator;
//use indicatif::ProgressIterator;
use glam::f64::DVec3;
use std::ops::Range;
use rand::Rng;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

fn main() -> io::Result<()>{
    //image
    let aspect_ratio: f64 = 16.0 / 9.0;
    let image_width: u32 = 1920;
    let camera: Camera = Camera::new(image_width, aspect_ratio);

    //World
    let world = HittableList{
        objects: vec![
            Box::new(Sphere{
                center: DVec3::new(0.0, -100.5, 1.0),
                radius: 100.0,
                material: Material::Diffuse {attenuation: DVec3::new(0.1, 0.8, 0.1)},
            }),

            Box::new(Sphere{
                center: DVec3::new(1.0, 0.0, 1.0),
                radius: 0.5,
                material: Material::Metal {attenuation: DVec3::new(0.95, 0.8, 0.15), fuzz: 0.0},
            }),
            Box::new(Sphere{
                center: DVec3::new(-1.5, -0.125, 1.75),
                radius: 0.25,
                material: Material::Diffuse {attenuation: DVec3::new(0.9, 0.1, 0.1)},
            }),
            Box::new(Sphere{
                center: DVec3::new(1.5, -0.125, 4.25),
                radius: 0.25,
                material: Material::Diffuse {attenuation: DVec3::new(0.9, 0.1, 0.1)},
            }),
            Box::new(Sphere{
                center: DVec3::new(2.0, -0.125, 4.25),
                radius: 0.25,
                material: Material::Diffuse {attenuation: DVec3::new(0.1, 0.1, 0.9)},
            }),
            Box::new(Sphere{
                center: DVec3::new(-0.2, -0.4375, 0.5),
                radius: 0.125,
                material: Material::Metal {attenuation: DVec3::new(0.1, 0.1, 0.8), fuzz: 0.5},
            }),
            Box::new(RectPrism{
                center: DVec3::new(-4.5, 2.0, 4.5),
                dimentions: DVec3::new(2., 0.2, 1.),
                material: Material::Metal {attenuation: DVec3::new(0.8, 0.5, 0.5), fuzz: 0.2},
            }),
            Box::new(RectPrism{
                center: DVec3::new(-4.5, 2.4, 4.5),
                dimentions: DVec3::new(2., 0.2, 1.),
                material: Material::Metal {attenuation: DVec3::new(0.8, 0.5, 0.5), fuzz: 0.2},
            }),
            Box::new(RectPrism{
                center: DVec3::new(-4.5, 2.8, 4.5),
                dimentions: DVec3::new(2., 0.2, 1.),
                material: Material::Metal {attenuation: DVec3::new(0.8, 0.5, 0.5), fuzz: 0.2},
            }),
            Box::new(RectPrism{
                center: DVec3::new(40.0, 40.0, 40.0),
                dimentions: DVec3::new(70.0, 100.0, 2.0),
                material: Material::Metal {attenuation: DVec3::new(0.6, 0.6, 0.6), fuzz: 0.0},
            }),
            Box::new(Sphere{
                center: DVec3::new(0.0, 0.0, 1.5),
                radius: 0.5,
                material: Material::Dielectric {attenuation: DVec3::new(0.98, 0.98, 0.98), index_of_refraction: 1.51},
            }),
            Box::new(RectPrism{
                center: DVec3::new(-3.5, 2.8, 7.0),
                dimentions: DVec3::new(1.0, 100.0, 1.0),
                material: Material::Metal {attenuation: DVec3::new(0.1, 0.1, 0.9), fuzz: 0.5},
            }),

            Box::new(RectPrism{
                center: DVec3::new(-1.0, -0.5, 1.0),
                dimentions: DVec3::new(0.4, 0.8, 0.4),
                material: Material::Diffuse {attenuation: DVec3::new(0.9, 0.1, 0.9)},
            }),

        ]
    };

    camera.render_to_disk(world)
}

struct Ray {
    origin: DVec3,
    direction: DVec3,
}

impl Ray {
    fn at(&self, t: f64) -> DVec3 {
        self.origin + t * self.direction
    }

    fn color<T>(&self, world: &T,  background_color: &DVec3, depth: u8) -> DVec3
    where T: Hittable,
    {

        let hit = world.hit(&self, &((0.001)..f64::INFINITY));

        if !hit.is_none(){
            let hit_record = hit.unwrap();
            let (new_ray, attenuation) = hit_record.material.scatter(&hit_record).unwrap();

            if depth > 0 {
                //Recursively return the final color
                return *attenuation*new_ray.color(world, background_color, depth-1);
            }
            return *attenuation;
        }
        let unit_direction: DVec3 = unit_vector(&self.direction);
        let a = 0.5 * (unit_direction.y + 1.0);
        return (1.0 - a) * DVec3::new(1.0, 1.0, 1.0) + a * DVec3::new(0.5, 0.7, 1.0);
    }
}

struct HitRecord {
    p: DVec3,
    normal: DVec3,
    incident_dir: DVec3,
    t: f64,
    material: Material,
    front_face: bool,
}

trait Hittable {
    fn hit(&self, r: &Ray, interval: &Range<f64>) -> Option<HitRecord>;
}

struct HittableList{
    objects: Vec<Box<dyn Hittable + Sync>>,
}


impl HittableList{
    fn clear(&mut self){
        self.objects = vec![]
    }

    fn add<T>(&mut self, object: T) where T: Hittable + Sync + 'static,{
        self.objects.push(Box::new(object));
    }
}

impl Hittable for HittableList{
    fn hit(&self, r: &Ray, interval: &Range<f64>) -> Option<HitRecord>{
        let mut closest_so_far: HitRecord = HitRecord{
            p: DVec3::ZERO,
            normal: DVec3::ZERO,
            incident_dir: DVec3::ZERO,
            t: 0.,
            material: Material::Diffuse {attenuation: DVec3::ZERO},
            front_face: true,
        };
        let mut hit_anything = false;

        for n in 0..(self.objects.len()) {
            //let hittable = *object;

            let hit = (*self.objects[n]).hit(r, interval);
            if !hit.is_none() {
                let record = hit.unwrap();
                if !hit_anything {
                    closest_so_far = record;
                    hit_anything = true;
                }
                else if record.t < closest_so_far.t{
                    closest_so_far = record;
                }
            }
        }
        if !hit_anything{
            return None;
        }
        return Some(closest_so_far);
    }
}

struct Sphere {
    center: DVec3,
    radius: f64,
    material: Material,
}

impl Hittable for Sphere {
    fn hit(&self, r: &Ray, interval: &Range<f64>) -> Option<HitRecord> {
        let oc: DVec3 = r.origin - self.center;

        let a: f64 = r.direction.dot(r.direction);
        let half_b: f64 = r.direction.dot(oc);
        let c: f64 = oc.dot(oc) - self.radius * self.radius;

        let discriminant = half_b * half_b - a * c;

        if discriminant < 0.0 {
            return None;
        }

        let sqrted: f64 = discriminant.sqrt();
        let p:DVec3;
        let t;

        let t1: f64 = (-half_b - sqrted) / a;
        if !interval.contains(&t1) {

            let t2: f64 = (-half_b + sqrted) / a;
            if !interval.contains(&t2) {
                return None;
            }
            p = r.at(t2);
            t=t2;

        }
        else{
            p = r.at(t1);
            t=t1;
        }
        let mut normal = (p - self.center) / self.radius;
        let front_face = r.direction.dot(normal) < 0.0;
        if !front_face {
            normal = -normal;
        }

        return Some(HitRecord{
            t,
            p,
            normal,
            incident_dir: r.direction,
            material: self.material.clone(),
            front_face,
        });
    }
}

struct RectPrism{
    center: DVec3,
    dimentions: DVec3,
    material: Material,
}

impl Hittable for RectPrism{
    fn hit(&self, r: &Ray, interval: &Range<f64>) -> Option<HitRecord>{

        let mut closest_so_far: HitRecord = HitRecord{
            p: DVec3::ZERO,
            normal: DVec3::ZERO,
            incident_dir: DVec3::ZERO,
            t: f64::MAX,
            material: Material::Diffuse {attenuation: DVec3::ZERO},
            front_face: true,
        };

        for x in 0..6 {
            //Check where the ray intersects plane defining each face
            let p;
            let normal;
            match x {
                0 => {p = self.center+DVec3::new(self.dimentions.x/2.0, 0.0, 0.0); normal = DVec3::X;}
                1 => {p = self.center+DVec3::new(0.0, self.dimentions.y/2.0, 0.0); normal = DVec3::Y;}
                2 => {p = self.center+DVec3::new(0.0, 0.0, self.dimentions.z/2.0); normal = DVec3::Z;}
                3 => {p = self.center-DVec3::new(self.dimentions.x/2.0, 0.0, 0.0); normal = DVec3::NEG_X;}
                4 => {p = self.center-DVec3::new(0.0, self.dimentions.y/2.0, 0.0); normal = DVec3::NEG_Y;}
                _ => {p = self.center-DVec3::new(0.0, 0.0, self.dimentions.z/2.0); normal = DVec3::NEG_Z;}
            }
            let t = intersects_plane(normal, p, r.direction, r.origin);

            let poi = r.at(t);

            //eliminate intersections outside of bounds of face
            let diff = self.center-self.dimentions/2.0;
            let sum = self.center+self.dimentions/2.0;

            let within_bounds: bool = match x%3{
                //yz
                0 => {poi.y >= diff.y && poi.y <= sum.y
                    && poi.z >= diff.z && poi.z <= sum.z}
                //xz
                1 => {poi.x >= diff.x && poi.x <= sum.x
                    && poi.z >= diff.z && poi.z <= sum.z}
                //xy
                _ => {poi.y >= diff.y && poi.y <= sum.y
                    && poi.x >= diff.x && poi.x <= sum.x }
            };
            if within_bounds {
                //eliminate intersections behind other intersections
                if t<closest_so_far.t && interval.contains(&t) {
                    let record = HitRecord{
                        p: poi,
                        normal,
                        incident_dir: r.direction,
                        t,
                        material: self.material.clone(),
                        front_face: r.direction.dot(normal) < 0.0,
                    };
                    closest_so_far = record;
                }

            }
        }
        //take closest one remaining or return none if it never hit
        if closest_so_far.t == f64::MAX{
            return None;
        }
        return Some(closest_so_far);
    }
}

fn intersects_plane(normal: DVec3, p: DVec3, m: DVec3, q: DVec3) -> f64{
    return (normal.dot(p)-normal.dot(q))/normal.dot(m);
}

struct Camera{
    image_width: u32,
    image_height: u32,
    max_value: u8,
    aspect_ratio: f64,
    center: DVec3,
    pixel_delta_u: DVec3,
    pixel_delta_v: DVec3,
    pixel00_loc: DVec3,
    samples_per_pixel: u8,
    background_color: DVec3,
    max_depth: u8,

}
impl Camera{
    fn new(image_width: u32, aspect_ratio: f64) -> Self{
        let max_value: u8 = 255;
        let image_height: u32 = ((image_width as f64) / aspect_ratio) as u32;
        if image_height < 1 {
            panic!("image height must be at least 1");
        }

        //Camera
        let focal_length: f64 = 1.0;
        let viewport_height: f64 = 2.0;
        let viewport_width = viewport_height * ((image_width as f64) / (image_height as f64));
        let camera_center: DVec3 = DVec3::ZERO;

        //Vectors along horizontal and down vertical edges of viewport
        let viewport_u: DVec3 = DVec3::new(viewport_width, 0.0, 0.0);
        let viewport_v: DVec3 = DVec3::new(0.0, -viewport_height, 0.0);

        //Horizontal and vertical delta vectors from pixel to pixel
        let pixel_delta_u: DVec3 = viewport_u / (image_width as f64);
        let pixel_delta_v: DVec3 = viewport_v / (image_height as f64);

        //Position of upper left pixel in viewport;
        let viewport_upper_left: DVec3 =
            DVec3::new(0.0, 0.0, focal_length) - viewport_u * 0.5 - viewport_v * 0.5;
        let pixel00_loc: DVec3 = viewport_upper_left + (pixel_delta_u + pixel_delta_v) * 0.5;

        return Camera{
            image_width,
            image_height,
            max_value,
            aspect_ratio,
            center: camera_center,
            pixel_delta_u,
            pixel_delta_v,
            pixel00_loc,
            samples_per_pixel: 64,
            background_color: DVec3::new(0.5, 0.7, 1.0),
            max_depth: 10,
        }

    }

    fn get_ray(&self, x:u32, y: u32) -> Ray{
        let mut rng = rand::thread_rng();
        let pixel_center = self.pixel00_loc
            + (x as f64) * self.pixel_delta_u
            + (y as f64) * self.pixel_delta_v;



        let (rand_x, rand_y) = (-0.5f64+rng.gen::<f64>(), -0.5f64+rng.gen::<f64>());

        let pixel_rand = (rand_x * self.pixel_delta_u)
            + (rand_y * self.pixel_delta_v);

        let pixel_center_rand = pixel_center + pixel_rand;

        let ray_direction = pixel_center_rand - self.center;
        return Ray { origin: self.center, direction: ray_direction };


    }

    fn render_to_disk<T>(&self, world: T) -> io::Result<()> where T: Hittable + std::marker::Sync,{
        //Render
        let mut file = File::create("img.ppm").expect("Could not create file");
        file.write_all(b"P3\n").expect("error writing to file");
        file.write_all(format!("{} {}\n", self.image_width, self.image_height).as_bytes()).expect(
            "error writing to file"
        );
        file.write_all(b"255\n").expect("error writing to file");



        let pixels = (0..self.image_height)
            .cartesian_product(0..self.image_width)
            .collect::<Vec<(u32, u32)>>()
            //.into_iter()
            .into_par_iter()
            .progress_count(
                self.image_height as u64
                    * self.image_width as u64,
            )
            .map(|(y, x)| {

                let scale_factor = (self.samples_per_pixel as f64).recip();

                let multisampled_pixel_color = (0..(self.samples_per_pixel)).into_iter().map(|_| {
                    self.get_ray(x, y)
                        .color(&world, &self.background_color, self.max_depth)
                    *scale_factor
                })
                    .sum::<DVec3>();
                //Apply linear to gamma transform
                let (r, g, b) = linear_to_gamma(multisampled_pixel_color.x, multisampled_pixel_color.y, multisampled_pixel_color.z);
                format!("{} {} {}", r as u8, g  as u8, b as u8)
            })
            .collect::<Vec<String>>()
            .join("\n");

        file.write_all(pixels.as_bytes()).expect("error writing to file");
        Ok(())
    }
}

fn unit_vector(vec: &DVec3) -> DVec3 {
    let len: f64 = (vec.x * vec.x + vec.y * vec.y + vec.z * vec.z).sqrt();
    *vec / len
}


fn rand_vec_on_hemisphere(normal: &DVec3) -> DVec3{
    let mut rng = rand::thread_rng();
    let mut rand_vec: DVec3;
    loop{
        //Generate random vec within the cube where -1 <= x, y, z <1
        rand_vec = DVec3::new(rng.gen::<f64>()*2.0-1.0, rng.gen::<f64>()*2.0-1.0, rng.gen::<f64>()*2.0-1.0);
        //Reject if outside the unit sphere
        let len = rand_vec.length_squared();
        if len<1.0{
            break;
        }
    }
    //Normalize
    let rand_vec = rand_vec.normalize();

    //Invert if in wrong hemisphere
    if rand_vec.dot(*normal)>=0.0 {
        return rand_vec;
    }
    return -1.*rand_vec;
}

#[derive(Clone)]
enum Material{
    Metal{attenuation: DVec3, fuzz: f64},
    Diffuse{attenuation: DVec3},
    Dielectric{attenuation: DVec3, index_of_refraction: f64},
}



impl Material{
    fn scatter(&self, hit_record: &HitRecord) -> Option<(Ray, &DVec3)>{
        match self {
            Material::Metal {attenuation, fuzz} =>{

                let direction = reflect(hit_record.incident_dir, hit_record.normal, fuzz);

                let new_ray = Ray{
                    origin: hit_record.p,
                    direction, };
                return Some((new_ray, attenuation));
            }

            Material::Diffuse {attenuation} =>{
                let mut direction = rand_vec_on_hemisphere(&hit_record.normal)+hit_record.normal;
                if is_near_zero(&direction) {
                    direction = hit_record.normal;
                }
                let new_ray = Ray{
                    origin: hit_record.p,
                    direction, };
                return Some((new_ray, attenuation));
            }

            Material::Dielectric {attenuation, index_of_refraction} => {
                let refraction_ratio: f64;
                if hit_record.front_face {    refraction_ratio = 1.0/index_of_refraction;  }
                else { refraction_ratio = *index_of_refraction; }

                let r_incoming = hit_record.incident_dir.normalize();
                let normal = hit_record.normal.normalize();

                let cos_theta = f64::min( (-r_incoming).dot(normal) , 1.0);
                let sin_theta = (1.0 - cos_theta*cos_theta).sqrt();

                let cannot_refract = refraction_ratio * sin_theta > 1.0;
                let direction;

                let mut rng = rand::thread_rng();

                if cannot_refract || reflectance(cos_theta, refraction_ratio) > rng.gen::<f64>() {
                    //reflect
                    direction = reflect(r_incoming, normal, &0.0);
                }else{
                    //can refract
                    direction = refract(r_incoming, normal, refraction_ratio);
                }

                //println!("Incoming: {}, \nNormal: {}, \nRefracted: {}",hit_record.incident_dir, hit_record.normal, refracted);
                //println!("ratio: {}, front-face: {}",refraction_ratio, hit_record.front_face);
                let scattered = Ray{
                    origin: hit_record.p,
                    direction,
                };
                Some((scattered, attenuation))
            }

        }
    }
}

fn reflectance(cosine: f64, ref_idx: f64) -> f64{
    //Schlick's approximation
    let r0: f64 = ((1.0-ref_idx) / (1.0+ref_idx)).powi(2);
    return r0 + (1.0-r0)*(1.0-cosine).powi(5);
}

fn linear_to_gamma(n1:f64, n2:f64, n3:f64) -> (f64, f64, f64){
    return (n1.sqrt()*255.0, n2.sqrt()*255.0, n3.sqrt()*255.0);
}

fn is_near_zero(vec3: &DVec3)->bool{
    let s = 1e-8;
    vec3.x < s && vec3.y < s && vec3.z < s
}

fn reflect(ray_direction: DVec3, normal: DVec3, fuzz: &f64) -> DVec3{
    if *fuzz > 1.0 || *fuzz <0.0{
        panic!("fuzz must be from 0 to 1")
    }
    let fuzz_vec = *fuzz*rand_vec_on_hemisphere(&normal);
    (ray_direction-2.0*ray_direction.dot(normal)*normal + fuzz_vec).normalize()
}

fn refract(r_incoming: DVec3, normal: DVec3, eta_over_eta_prime: f64) -> DVec3{
    let cos_theta = f64::min( (-r_incoming).dot(normal) , 1.0);
    //find the components of r prime which are parallel and perpendicular to normal
    let r_prime_perpendicular = eta_over_eta_prime*(r_incoming + cos_theta*normal);
    let r_prime_parallel = -(1.0-r_prime_perpendicular.length_squared()).abs().sqrt()*normal;

    return r_prime_perpendicular + r_prime_parallel;
}