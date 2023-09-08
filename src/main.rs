use std::fs::File;
use std::io::prelude::*;
use itertools::Itertools;
use indicatif::ProgressIterator;
use glam::f64::DVec3;

fn main() {
    //image
    let aspect_ratio: f64 = 16.0 / 9.0;
    let image_width: u32 = 200;
    let image_height: u32 = ((image_width as f64) / aspect_ratio) as u32;
    if image_height < 1 {
        panic!("iumage height must be at least 1");
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
    let pixel00_location: DVec3 = viewport_upper_left + (pixel_delta_u + pixel_delta_v) * 0.5;

    //Render
    let mut file = File::create("img.ppm").expect("Could not create file");
    file.write_all(b"P3\n").expect("error writing to file");
    file.write_all(format!("{} {}\n", image_width, image_height).as_bytes()).expect(
        "error writing to file"
    );
    file.write_all(b"255\n").expect("error writing to file");

    let pixels = (0..image_height)
        .cartesian_product(0..image_width)
        .progress_count((image_width as u64) * (image_height as u64))
        .map(|(y, x)| {
            let pixel_center =
                pixel00_location + (x as f64) * pixel_delta_u + (y as f64) * pixel_delta_v;
            let ray_direction = pixel_center - camera_center;
            let ray: Ray = Ray { origin: camera_center, direction: ray_direction };

            let pixel_color_float = ray_color(&ray);
            let pixel_color: (u8, u8, u8) = (
                (pixel_color_float.x * 255.0) as u8,
                (pixel_color_float.y * 255.0) as u8,
                (pixel_color_float.z * 255.0) as u8,
            );
            format!("{} {} {}", pixel_color.0, pixel_color.1, pixel_color.2)
        })
        .join("\n");

    file.write_all(pixels.as_bytes()).expect("error writing to file");
}

struct Ray {
    origin: DVec3,
    direction: DVec3,
}

impl Ray {
    fn at(&self, t: f64) -> DVec3 {
        self.origin + t * self.direction
    }
}

struct HitRecord {
    p: DVec3,
    normal: DVec3,
    t: f64,
}

trait Hittable {
    fn hit(&self, r: &Ray, ray_t_min: f64, ray_t_max: f64, rec: &mut HitRecord) -> bool;
}

struct Sphere {
    center: DVec3,
    radius: f64,
}

impl Hittable for Sphere {
    fn hit(&self, r: &Ray, ray_t_min: f64, ray_t_max: f64, rec: &mut HitRecord) -> bool {
        let oc: DVec3 = r.origin - self.center;

        let a: f64 = r.direction.clone().dot(r.direction.clone());
        let half_b: f64 = r.direction.clone().dot(oc.clone());
        let c: f64 = oc.clone().dot(oc.clone()) - self.radius * self.radius;

        let discriminant = half_b * half_b - a * c;

        if discriminant < 0.0 {
            return false;
        }

        let sqrted: f64 = discriminant.sqrt();

        let t1: f64 = (-half_b - sqrted) / a;
        if t1 <= ray_t_min || t1 >= ray_t_max {
            let t2: f64 = (-half_b + sqrted) / a;
            if t2 <= ray_t_min || t2 >= ray_t_max {
                return false;
            }
            rec.t = t2;
            rec.p = r.at(rec.t);
            rec.normal = (rec.p - self.center) / self.radius;
            return true;
        }
        rec.t = t1;
        rec.p = r.at(rec.t);
        rec.normal = (rec.p - self.center) / self.radius;
        return true;
    }
}

fn unit_vector(vec: &DVec3) -> DVec3 {
    let len: f64 = (vec.x * vec.x + vec.y * vec.y + vec.z * vec.z).sqrt();
    vec.clone() / len
}

fn hit_sphere(center: &DVec3, radius: f64, r: &Ray) -> f64 {
    let oc: DVec3 = r.origin - *center;

    let a: f64 = r.direction.clone().dot(r.direction.clone());
    let half_b: f64 = r.direction.clone().dot(oc.clone());
    let c: f64 = oc.clone().dot(oc.clone()) - radius * radius;

    let discriminant = half_b * half_b - a * c;

    if discriminant < 0.0 {
        -1.0
    } else {
        (-half_b - discriminant.sqrt()) / a
    }
}

fn ray_color(ray: &Ray) -> DVec3 {
    let center: DVec3 = DVec3::new(0.0, 0.0, 1.0);
    let radius: f64 = 0.5;

    let t = hit_sphere(&center, radius, ray);

    if t > 0.0 {
        let n: DVec3 = (ray.at(t) - center) / radius;
        //println!("n: {}",n);
        let color: DVec3 = 0.5 * DVec3::new(n.x + 1.0, n.y + 1.0, n.z + 1.0);
        //println!("color: {}",color);
        return color;
    }

    let unit_direction: DVec3 = unit_vector(&ray.direction);
    let a = 0.5 * (unit_direction.y + 1.0);
    (1.0 - a) * DVec3::new(1.0, 1.0, 1.0) + a * DVec3::new(0.5, 0.7, 1.0)
}
