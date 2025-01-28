# Ray Tracer in Rust

Welcome to **Ray Tracer in Rust**, a project inspired by the "Ray Tracing in One Weekend" book series ([official site](https://raytracing.github.io/books/RayTracingInOneWeekend.html)) and Chris Biscardi's "Ray Tracing in Rust" YouTube video. This project explores the fundamentals of computer graphics, focusing on rendering realistic scenes with light and geometry.

## Project Overview

This **basic ray tracer** is written in **Rust**, leveraging the language's performance and safety features to create visually compelling renderings. It processes a list of "hittable" objects, simulating light interactions such as reflections and shadows, and outputs the final render as a `.ppm` image file. 

![Rendered Image](https://github.com/TomChiavegato/ray_tracer/assets/129907786/340fe0b1-90c1-4a8e-9b80-0b475e7ecdd5)

### Key Features
- **Custom Object-Oriented Design**: Implements a modular approach with "hittable" objects that allow extensibility for additional shapes and materials.
- **High-Quality Rendering**: Supports features like anti-aliasing, reflections, and multiple light bounces for realistic output.
- **Optimized Performance**: Leverages Rust's efficient memory management for fast rendering and smooth execution.

### Learning Goals
This project allowed me to:
- Deepen my understanding of **ray tracing algorithms** and graphics programming.
- Explore **Rust programming**, including concepts like lifetimes, ownership, and crates for math and image handling.
- Develop a hands-on appreciation for **light physics** and its interaction with objects in 3D space.

## Getting Started

### Prerequisites
- **Rust**: Ensure you have Rust installed. [Install Rust](https://www.rust-lang.org/tools/install)
- A code editor or IDE with Rust support (e.g., Visual Studio Code or IntelliJ IDEA).

### Running the Project
1. Clone the repository:
   ```bash
   git clone https://github.com/TomChiavegato/ray_tracer.git
   cd ray_tracer
   ```
2. Build and run the project:
   ```bash
   cargo run
   ```
3. Find the generated `.ppm` file in the project directory and open it using an image viewer that supports `.ppm` (e.g., GIMP).

## Features to Implement (Roadmap)
- **Code Refactor**: Organize the codebase into modular files and a `lib.rs` for better maintainability.
- **Defocus Blur (Depth of Field)**: Add support for realistic blurring effects.

## Lessons and Personal Contributions
While inspired by existing resources, I **independently wrote the majority of the code** and customized the implementation to my own understanding. This included:
- Experimenting with different **sampling methods** for anti-aliasing.
- Structuring the project to adhere to Rust's idiomatic practices.
- Implementing rectangular prisim objects.

## Why This Project Stands Out
This ray tracer highlights:
- My ability to **learn independently** by tackling complex concepts like light physics and 3D rendering.
- My capability to **adapt and extend** existing methodologies into personalized solutions.
- My technical strengths in **Rust programming**, a language widely recognized for its robustness and performance.

## Contact
If you have feedback or would like to discuss this project further, feel free to connect with me on [LinkedIn](https://linkedin.com/in/tom-chiavegato) or reach out via email: tomchiavegato@gmail.com.
