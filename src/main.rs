use std::fs::File;
use std::io::prelude::*;
use itertools::Itertools;
use indicatif::ProgressIterator;
use glam::f64::DVec3;

fn main() {
    //image
    let aspect_ratio: f64 = 16.0/9.0;
    let image_width: u32 = 800;
    let image_height: u32 = (image_width as f64/aspect_ratio) as u32;
    if image_height <1{
        panic!("iumage height must be at least 1");
    }

    //Camera
    let focal_length: f64 = 1.0;
    let viewport_height: f64 = 2.0;
    let viewport_width = viewport_height * (image_width as f64/image_height as f64);
    let camera_center: DVec3 = DVec3::new(0., 0., 0.);

    //Vectors along horizontal and down vertical edges of viewport
    let viewport_u: DVec3 = DVec3::new(viewport_width, 0., 0.);
    let viewport_v: DVec3 = DVec3::new(0., -viewport_height, 0.);

    //Horizontal and vertical delta vectors from pixel to pixel
    let pixel_delta_u: DVec3 = viewport_u / image_width as f64;
    let pixel_delta_v: DVec3 = viewport_v / image_height as f64;

    //Position of upper left pixel in viewport;
    let viewport_upper_left: DVec3 = DVec3::new(0., 0., focal_length) - viewport_u*0.5 - viewport_v*0.5;
    let pixel00_location: DVec3 = viewport_upper_left + (pixel_delta_u+ pixel_delta_v)*0.5;


    //Render
    let mut file = File::create("img.ppm")
        .expect("Could not create file");
    file.write_all(b"P3\n")
        .expect("error writing to file");
    file.write_all(format!("{} {}\n", image_width, image_height).as_bytes())
        .expect("error writing to file");
    file.write_all(b"255\n")
        .expect("error writing to file");

    let pixels = (0..image_height)
        .cartesian_product(0..image_width)
        .progress_count(image_width as u64 * image_height as u64)
        .map(
            |(y, x)| {
                let pixel_center = pixel00_location + x as f64*pixel_delta_u + y as f64*pixel_delta_v;
                let ray_direction = pixel_center - camera_center;
                let ray: Ray = Ray{ origin: camera_center, direction: ray_direction,};

                let pixel_color_float = ray_color(&ray);
                let pixel_color: (u8,u8,u8) = ((pixel_color_float.x*255.0) as u8, (pixel_color_float.y*255.0) as u8, (pixel_color_float.z*255.0) as u8);
                format!("{} {} {}", pixel_color.0, pixel_color.1, pixel_color.2)
        
        })
        .join("\n");
    
    file.write_all(pixels.as_bytes())
    .expect("error writing to file");
}



struct Ray{
    origin: DVec3,
    direction: DVec3,
}

impl Ray{
    fn at(&self, t: f64) -> DVec3 {
        self.origin+t*self.direction
    }
}

fn unit_vector(vec: &DVec3)->DVec3{
    let len: f64 = (vec.x*vec.x+vec.y*vec.y+vec.z*vec.z).sqrt();
    vec.clone()/len
}

fn hit_sphere(center: &DVec3, radius: f64, r: &Ray) -> bool{
    
    let oc: DVec3 = r.origin-*center;

    let a: f64 = r.direction.clone().dot(r.direction.clone());
    let b: f64 = 2. * r.direction.clone().dot(oc.clone());
    let c: f64 = oc.clone().dot(oc.clone()) - radius*radius;

    let discriminant = b*b - 4. * a * c;
    
    discriminant >= 0.
}

fn ray_color (ray: &Ray)->DVec3{
    if  hit_sphere(&DVec3::new(0.,0.,1.), 0.5, ray) {
        return DVec3::new(1. , 0. ,0. );
    }
    let unit_direction:DVec3 = unit_vector(&ray.direction);
    let a = 0.5*(unit_direction.y + 1.0);
    (1.0-a)*DVec3::new(1.0,1.0,1.0)+a*DVec3::new(0.5,0.7,1.0)
}