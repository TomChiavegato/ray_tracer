# ray_tracer
<p>
  Basic ray tracer made in Rust by following the c++ "ray tracing in a weekend" book https://raytracing.github.io/books/RayTracingInOneWeekend.html#surfacenormalsandmultipleobjects.
I got the idea from Chris Biscardi's YouTube video "raytracing in rust" where he did the same thing.I used used some of the same packages and project structure he talked about but coded
most of the project without looking at the video.
The ray tracer takes a list of "hittable" objects and then writes the render to a .ppm image file.
It currently all runs on one thread so is pretty slow but I plan to change this in the future.
</p>

![img](https://github.com/TomChiavegato/ray_tracer/assets/129907786/2e3ce5e7-634f-488c-af68-db292eb9c7d1)

<h3>Still need to do:</h4>
<ul>
  <li>Organize code into seperate files and refactor into a lib.rs</li>
  <li>Add parallelism so </li>
  <li>Add more shapes</li>
</ul>



