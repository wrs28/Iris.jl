using LinearAlgebra

## 1D

p1 = Point(3)
p2 = Point{Cartesian}(3)
p1 == p2

p1*.3333
.3333*p1
p1/3
3\p1

p1+p2
p1-p2

norm(p1)

ndims(p1)
p1[1]
p1<4
p1.x


## 2D

p1 = Point(3,4)
p2 = Point{Cartesian}(3,4)
p1 == p2

p1*.3333
.3333*p1
p1/3
3\p1

p1+p2
p1-p2

norm(p1)

ndims(p1)
p1[1]

p1.x
p1.y
p1.r
p1.ϕ

p1 = Point{Polar}(3,4)
p2 = Point{Polar}(4,5)

.3333*p1
p1/3
3\p1

p1+p2
p1-p2

norm(p1)

ndims(p1)
p1[1]

p1.x
p1.y
p1.r
p1.ϕ
p1.s

p1+p2


## 3D

p1 = Point{Polar}(3,4,5)
p2 = Point{Cartesian}(3,4,6)
p3 = Point{Spherical}(3,4,6)

p1*.3333
.3333*p1
p1/3
3\p1

p1+p2
p1+p3
p1-p2
p1-p3

norm(p1)
norm(p2)
norm(p3)

ndims(p1)
p1[1]

p1.x
p1.y
p1.z
p1.r
p1.ϕ
p1.θ

p1 = Point{Spherical}(3,4,6)
p2 = Point{Spherical}(4,5,6)

.3333*p1
p1/3
3\p1

p1+p2
p1-p2

norm(p1)

ndims(p1)
p1[1]

p1.x
p1.y
p1.z
p1.r
p1.ϕ
p1.s
p1.θ

p1+p2
