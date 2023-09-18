use std::fs::File;
use std::io;
use std::io::prelude::*;
use itertools::Itertools;
use indicatif::ProgressIterator;
use glam::f64::DVec3;
use std::ops::Range;
use rand::Rng;

fn main() -> io::Result<()>{
    let rand = rand_vec_on_hemisphere(&DVec3::new(1.,1.,1.));
    println!("{}",rand);
    //image
    let aspect_ratio: f64 = 16.0 / 9.0;
    let image_width: u32 = 600;

    let camera: Camera = Camera::new(image_width, aspect_ratio);

    //World
    let world = HittableList{
        objects: vec![
            Box::new(Sphere{
                center: DVec3::new(0.0, 0.0, 1.0),
                radius: 0.5,
                color: DVec3::new(0.5, 0.5, 0.5),
            }),
            Box::new(Sphere{
                center: DVec3::new(0.0, -100.5, 1.0),
                radius: 100.0,
                color: DVec3::new(0.5, 0.5, 0.5),
            }),
        ]
    };

    camera.render_to_disk(&world);
    Ok(())
}

struct Ray {
    origin: DVec3,
    direction: DVec3,
}

impl Ray {
    fn at(&self, t: f64) -> DVec3 {
        self.origin + t * self.direction
    }

    fn color<T>(&self, world: &T, color: DVec3, background_color: &DVec3, depth: u8) -> DVec3
    where T: Hittable,
    {
        let hit = world.hit(&self, &((0.)..f64::INFINITY));

        if !hit.is_none(){
            let hit_record = hit.unwrap();
            let new_ray = Ray{
                origin: hit_record.p,
                direction: rand_vec_on_hemisphere(&hit_record.normal), };
            if depth > 0 {
                return new_ray.color(world, color*hit_record.color, background_color, depth-1);
            }
            return hit_record.color*hit_record.color
        }
        let unit_direction: DVec3 = unit_vector(&self.direction);
        let a = 0.5 * (unit_direction.y + 1.0);
        return (1.0 - a) * DVec3::new(1.0, 1.0, 1.0) + a * DVec3::new(0.5, 0.7, 1.0);
    }
}

struct HitRecord {
    p: DVec3,
    normal: DVec3,
    t: f64,
    color: DVec3,
}

trait Hittable {
    fn hit(&self, r: &Ray, interval: &Range<f64>) -> Option<HitRecord>;
}

struct HittableList{
    objects: Vec<Box<dyn Hittable>>,
}
impl HittableList{
    fn clear(&mut self){
        self.objects = vec![]
    }

    fn add<T>(&mut self, object: T) where T: Hittable + 'static,{
        self.objects.push(Box::new(object));
    }
}

impl Hittable for HittableList{
    fn hit(&self, r: &Ray, interval: &Range<f64>) -> Option<HitRecord>{
        let mut closest_so_far: HitRecord = HitRecord{
            p: DVec3::new(0., 0., 0.),
            normal: DVec3::new(0., 0., 0.),
            t: 0.,
            color: DVec3::new(0., 0., 0.),
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
    color: DVec3,
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

        let t1: f64 = (-half_b - sqrted) / a;
        if !interval.contains(&t1) {
            let t2: f64 = (-half_b + sqrted) / a;
            if !interval.contains(&t2) {
                return None;
            }
            let p = r.at(t2);
            return Some(HitRecord{
                t: t2,
                p,
                normal: (p - self.center) / self.radius,
                color: self.color,
            });
        }
        let p = r.at(t1);
        return Some(HitRecord{
            t: t1,
            p,
            normal: (p - self.center) / self.radius,
            color: self.color,
        });
    }
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
    depth: u8,

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
        let camera_center: DVec3 = DVec3::new(0.0, 0.0, 0.0);

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
            samples_per_pixel: 16,
            background_color: DVec3::new(0.5, 0.7, 1.0),
            depth: 30,
        }

    }

    fn get_ray(&self, x:u32, y: u32) -> Ray{
        let mut rng = rand::thread_rng();
        let pixel_center = self.pixel00_loc
            + (x as f64) * self.pixel_delta_u
            + (y as f64) * self.pixel_delta_v;



        let (randX, randY) = (-0.5f64+rng.gen::<f64>(), -0.5f64+rng.gen::<f64>());

        let pixel_rand = (randX * self.pixel_delta_u)
            + (randY * self.pixel_delta_v);

        let pixel_center_rand = pixel_center + pixel_rand;

        let ray_direction = pixel_center_rand - self.center;
        return Ray { origin: self.center, direction: ray_direction };


    }

    fn render_to_disk<T>(&self, world: &T) -> io::Result<()> where T: Hittable,{
        //Render
        let mut file = File::create("img.ppm").expect("Could not create file");
        file.write_all(b"P3\n").expect("error writing to file");
        file.write_all(format!("{} {}\n", self.image_width, self.image_height).as_bytes()).expect(
            "error writing to file"
        );
        file.write_all(b"255\n").expect("error writing to file");



        let pixels = (0..self.image_height)
            .cartesian_product(0..self.image_width)
            .progress_count((self.image_width as u64) * (self.image_height as u64))
            .map(|(y, x)| {

                let scale_factor = (self.samples_per_pixel as f64).recip();

                let multisampled_pixel_color = (0..(self.samples_per_pixel)).map(|_| {
                    self.get_ray(x, y)
                        .color(world, DVec3::new(1.0, 1.0, 1.0), &self.background_color, self.depth)
                        *255.
                    *scale_factor
                })
                    .sum::<DVec3>();


                format!("{} {} {}", multisampled_pixel_color.x as u8, multisampled_pixel_color.y  as u8, multisampled_pixel_color.z as u8)
            })
            .join("\n");

        file.write_all(pixels.as_bytes()).expect("error writing to file");
        Ok(())
    }
}

fn unit_vector(vec: &DVec3) -> DVec3 {
    let len: f64 = (vec.x * vec.x + vec.y * vec.y + vec.z * vec.z).sqrt();
    *vec / len
}

fn hit_sphere(center: &DVec3, radius: f64, r: &Ray) -> f64 {
    let oc: DVec3 = r.origin - *center;

    let a: f64 = r.direction.dot(r.direction);
    let half_b: f64 = r.direction.dot(oc);
    let c: f64 = oc.dot(oc) - radius * radius;

    let discriminant = half_b * half_b - a * c;

    if discriminant < 0.0 {
        -1.0
    } else {
        (-half_b - discriminant.sqrt()) / a
    }
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