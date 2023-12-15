use nalgebra::base::*;

#[derive(Copy, Clone)]
pub struct Circle {
    pub center: Vector2<f32>,
    pub radius: f32,
}

#[derive(Copy, Clone)]
pub struct Line {
    pub normal: Vector2<f32>,
    pub factor: f32,
}

pub enum ExtendedCircle {
    Circle(Circle),
    Line(Line),
}

impl ExtendedCircle {
    pub fn is_vertex_in(&self, vert: &Vector2<f32>) -> bool {
        match self {
            ExtendedCircle::Circle(circle) => circle.is_vertex_in(vert),
            ExtendedCircle::Line(line) => line.is_vertex_in(vert),
        }
    }
}

impl Circle {
    pub fn new(center: Vector2<f32>, radius: f32) -> Circle {
        Circle { center, radius }
    }
    pub fn is_vertex_in(&self, vert: &Vector2<f32>) -> bool {
        (self.center - vert).norm() - self.radius <= 0.
    }
}

impl Line {
    pub fn new(normal: Vector2<f32>, factor: f32) -> Line {
        Line { normal, factor }
    }
    pub fn is_vertex_in(&self, vert: &Vector2<f32>) -> bool {
        self.normal.dot(&vert) - self.factor <= 0.
    }
}
