-- lua script applying available C++ functions

-- vec2d
print("\nvec2d related:\n")

local v0 = vec2d.new()
local v1 = vec2d.new(1,1)
local v2 = vec2d.new(1.5,-2)
local v3 = vec2d.new(-v2)
local v4 = vec2d.new(2.5,-1)

print("v0:",v0, "empty ctor")
print("v1:",v1, "component ctor")
print("v2:",v2, "component ctor")
print("v3:",v3, "vector ctor v3 = -v2", "unary minus")
print("v4:",v4, "component ctor v4 = v1 + v3 (by hand)")

assert(v4 == v1 + v2)
assert(v1 == v4 - v2)
if (v4 == v1 + v2 and v1 == v4 - v2) then
    print ("v4 == v1 + v2", "comparison, addition and subtraction")
end

v5 = 3 * v4
assert(3 * v4 == v4 * 3)
assert(v5 / 3 == v4)
print("3 * v4:", 3 * v4, "scalar multiplication")
print("v5 / 3:", v5 / 3, "scalar division")

print("dot(v1, v2):", dot(v1, v2), "dot product")

print("sq_nrm(v1): ", sq_nrm(v1), "squared norm")
print("nrm(v1): ", nrm(v1), "norm")
assert(math.abs(sq_nrm(v1)-2) < eps)
assert(math.abs(nrm(v1)-math.sqrt(2)) < eps)

print("unitized(v1): ", unitized(v1), "unitized")
print("math.abs( nrm(unitized(v1))-1.0 ):", math.abs(nrm(unitized(v1))-1.0))
assert(math.abs(nrm(unitized(v1))-1.0) < eps)

print("inv(v1):", inv(v1), "inverse")
print("math.abs(nrm(v1)*nrm(inv(v1))-1):", math.abs(nrm(v1)*nrm(inv(v1))-1))
assert(math.abs(nrm(v1)*nrm(inv(v1))-1) < eps)

print("wdg(v1,v2):", wdg(v1,v2))

local ps1 = pscalar2d.new()
local ps2 = pscalar2d.new(3.7)
assert(ps1 == pscalar2d.new())
print("ps1:", ps1)
print("ps2:", ps2)

local phi = angle(v1,v2)
print("phi (rad): ", phi)
print("phi (deg): ", rad_to_deg(phi))
assert(math.abs(angle(vec2d.new(1,0),vec2d.new(1,1))-math.pi/4) < eps)

-- mvec2d_e
print("\nmvec2d_e related:\n")

local mve0 = mvec2d_e.new()
local mve1 = mvec2d_e.new(1,1)
local mve2 = mvec2d_e.new(1.5,-2)
local mve3 = mvec2d_e.new(-mve2)
local mve4 = mvec2d_e.new(2.5,-1)
local mve5 = mvec2d_e.new(scalar2d.new(-5))
local mve6 = mvec2d_e.new(pscalar2d.new(-6))
local mve7 = mvec2d_e.new(scalar2d.new(-7),pscalar2d.new(-8))

print("mve0:",mve0, "empty ctor")
print("mve1:",mve1, "component ctor")
print("mve2:",mve2, "component ctor")
print("mve3:",mve3, "even multivector ctor v3 = -v2", "unary minus")
print("mve4:",mve4, "component ctor v4 = v1 + v3 (by hand)")
print("mve5:",mve5, "ctor mvec2d(scalar2d(-5))")
print("mve6:",mve6, "ctor mvec2d(pscalar2d(-6))")
print("mve7:",mve7, "ctor mvec2d(scalar2d(-7),pscalar2d(-8))")

assert(mve4 == mve1 + mve2)
assert(mve1 == mve4 - mve2)
if (mve4 == mve1 + mve2 and mve1 == mve4 - mve2) then
    print ("mve4 == mve1 + mve2", "comparison, addition and subtraction")
end

mve5 = 3 * mve4
assert(3 * mve4 == mve4 * 3)
assert(mve5 / 3 == mve4)
print("mv5 = 3 * mve4:", 3 * mve4, "scalar multiplication")
print("mve5 / 3:", mve5 / 3, "scalar division")

print("gr0(mve2)", gr0(mve2))
print("gr2(mve2)", gr2(mve2))
assert(gr0(mve2) == scalar2d.new(1.5))
assert(gr2(mve2) == pscalar2d.new(-2))

print("angle_to_re(mve1) - (rad, deg):", angle_to_re(mve1), rad_to_deg(angle_to_re(mve1)))
print("angle_to_re(mve2) - (rad, deg):", angle_to_re(mve2), rad_to_deg(angle_to_re(mve2)))

-- inverses
print("\ninverses:\n")
v = vec2d.new(13,5)
vi = inv(v)
print("v:", v)
print("vi:", vi)
print("v*vi:", v*vi)
assert(math.abs(to_val(gr0(v*vi)) - 1) < eps)
assert(math.abs(to_val(gr2(v*vi)) - 0) < eps)
print()

vc = mvec2d_e.new(13,5)
vci = inv(vc)
print("vc:", vc)
print("vci:", vci)
print("vc*vci:", vc*vci)
assert(math.abs(to_val(gr0(vc*vci)) - 1) < eps)
assert(math.abs(to_val(gr2(vc*vci)) - 0) < eps)
print()

vm = mvec2d.new(13,-27,3,5)
vmi = inv(vm)
print("vm:", vm)
print("vmi:", vmi)
print("vm*vmi:", vm*vmi)
assert(math.abs(to_val(gr0(vm*vmi)) - 1) < eps)
assert(nrm(gr1(vm*vmi)) < eps)
assert(math.abs(to_val(gr2(vm*vmi)) - 0) < eps)
print()

-- mvec2d
print("\nmvec2d related:\n")

local mv0 = mvec2d.new()
local mv1 = mvec2d.new(1,1,1,1)
local mv2 = mvec2d.new(1.5,-2,1,2)
local mv3 = mvec2d.new(-mv2)
local mv4 = mvec2d.new(2.5,-1,2,3)
local mv5 = mvec2d.new(scalar2d.new(-5))
local mv6 = mvec2d.new(pscalar2d.new(-6))
local mv7 = mvec2d.new(scalar2d.new(-7),pscalar2d.new(-8))
local mv8 = mvec2d.new(vec2d.new(-2,1))
local mv9 = mvec2d.new(mvec2d_e.new(2.5,3))

print("mv0:",mv0, "empty ctor")
print("mv1:",mv1, "component ctor")
print("mv2:",mv2, "component ctor")
print("mv3:",mv3, "multivector ctor v3 = -v2", "unary minus")
print("mv4:",mv4, "component ctor v4 = v1 + v3 (by hand)")
print("mv5:",mv5, "ctor mvec2d(scalar2d(-5))")
print("mv6:",mv6, "ctor mvec2d(pscalar2d(-6))")
print("mv7:",mv7, "ctor mvec2d(scalar2d(-7),pscalar2d(-8))")
print("mv8:",mv8, "ctor mvec2d(vec2d(-2,1))")
print("mv9:",mv9, "ctor mvec2d(mvec2d_e(2.5,3))")

assert(mv4 == mv1 + mv2)
assert(mv1 == mv4 - mv2)
if (mv4 == mv1 + mv2 and mv1 == mv4 - mv2) then
    print ("mv4 == mv1 + mv2", "comparison, addition and subtraction")
end

mv5 = 3 * mv4
assert(3 * mv4 == mv4 * 3)
assert(mv5 / 3 == mv4)
print("mv5 = 3 * mv4:", 3 * mv4, "scalar multiplication")
print("mv5 / 3:", mv5 / 3, "scalar division")

print("gr0(mv2)", gr0(mv2))
print("gr1(mv2)", gr1(mv2))
print("gr2(mv2)", gr2(mv2))
assert(gr0(mv2) == scalar2d.new(1.5))
assert(gr1(mv2) == vec2d.new(-2,1))
assert(gr2(mv2) == pscalar2d.new(2))

print("sq_nrm(mv1), nrm(mv1)", sq_nrm(mv1), nrm(mv1))
print("rev(mv1), conj(mv1), unitized(mv1)", rev(mv1), conj(mv1), unitized(mv1))

print("inv(mv2), mv2*inv(mv2)", inv(mv2), mv2*inv(mv2))

-- handling of pscalar2d
print("I_2d*I_2d:", I_2d*I_2d, "sq_nrm(I_2d):", sq_nrm(I_2d), "nrm(I_2d):",nrm(I_2d))
I = pscalar2d.new(3)
print("I:", I, "I*inv(I):", I*inv(I))
assert(I*inv(I) == 1)

print("I+I:", I+I,"I-I:", I-I,"I/3:", I/3,"to_val(I):", to_val(I))

-- 2d_ops

-- projections, rejections and reflections
print("project_onto(vec2d.new(3,7),vec2d.new(1,1))", project_onto(vec2d.new(3,7),vec2d.new(1,1)))
print("reject_from(vec2d.new(3,7),vec2d.new(1,1))", reject_from(vec2d.new(3,7),vec2d.new(1,1)))

print("reflect_on_vec(vec2d.new(-2,2),vec2d.new(0,1))", reflect_on_vec(vec2d.new(-2,2),vec2d.new(0,1)))
print("reflect_on_hyp(vec2d.new(-2,2),vec2d.new(0,1))", reflect_on_hyp(vec2d.new(-2,2),vec2d.new(0,1)))

-- geometric products
print("mvec2d.new(2,-3,5,7)*mvec2d.new(1,2,-5,12):",mvec2d.new(2,-3,5,7)*mvec2d.new(1,2,-5,12))
print()
print("dot(vec2d.new(2,-3),vec2d.new(1,2)):",dot(vec2d.new(2,-3),vec2d.new(1,2)))
print("wdg(vec2d.new(2,-3),vec2d.new(1,2)):",to_val(wdg(vec2d.new(2,-3),vec2d.new(1,2))))
print("angle(vec2d.new(2,-3),vec2d.new(1,2)) [deg]:",rad_to_deg(angle(vec2d.new(2,-3),vec2d.new(1,2))))
print("nrm(vec2d.new(2,-3)):",nrm(vec2d.new(2,-3)))
print("nrm(vec2d.new(1,2)):",nrm(vec2d.new(1,2)))
print("nrm(vec2d.new(2,-3))*nrm(vec2d.new(1,2))", nrm(vec2d.new(2,-3))*nrm(vec2d.new(1,2)))
print("sin(angle)*nrm(vec2d.new(2,-3))*nrm(vec2d.new(1,2)))",
math.sin(angle(vec2d.new(2,-3),vec2d.new(1,2)))*nrm(vec2d.new(2,-3))*nrm(vec2d.new(1,2)))
print()
print("vec2d.new(2,-3)*vec2d.new(1,2):",vec2d.new(2,-3)*vec2d.new(1,2))
print()
print("I_2d*mvec2d.new(2,-3,5,7)",I_2d*mvec2d.new(2,-3,5,7))
print("dual(mvec2d.new(2,-3,5,7))",dual(mvec2d.new(2,-3,5,7)))
print()
print("I_2d*mvec2d_e.new(2,7)",I_2d*mvec2d_e.new(2,7))
print("dual(mvec2d_e.new(2,7))",dual(mvec2d_e.new(2,7)))
print()
print("I_2d*vec2d.new(2,7)",I_2d*vec2d.new(2,7))
print("dual(vec2d.new(2,7))",dual(vec2d.new(2,7)))
print()
print("mvec2d.new(2,-3,5,7)*I_2d",mvec2d.new(2,-3,5,7)*I_2d)
print("dual(mvec2d.new(2,-3,5,7)*I_2d)",dual(mvec2d.new(2,-3,5,7)))
print()
print("mvec2d_e.new(2,7)*I_2d",mvec2d_e.new(2,7)*I_2d)
print("dual(mvec2d_e.new(2,7))",dual(mvec2d_e.new(2,7)))
print()
print("vec2d.new(2,7)*I_2d",vec2d.new(2,7)*I_2d)
print("dual(vec2d.new(2,7))",dual(vec2d.new(2,7)))
print()
print("dual(scalar2d.new(5))", dual(scalar2d.new(5)))
print("dual(pscalar2d.new(3))", dual(pscalar2d.new(3)))
print()
print("vec2d.new(2,7)*mvec2d.new(1,2,-5,12):",vec2d.new(2,7)*mvec2d.new(1,2,-5,12))
print("mvec2d.new(0,2,7,0)*mvec2d.new(1,2,-5,12):",mvec2d.new(0,2,7,0)*mvec2d.new(1,2,-5,12))
print()
print("mvec2d.new(1,2,-5,12)*vec2d.new(2,7):",mvec2d.new(1,2,-5,12)*vec2d.new(2,7))
print("mvec2d.new(1,2,-5,12)*mvec2d.new(0,2,7,0):",mvec2d.new(1,2,-5,12)*mvec2d.new(0,2,7,0))
print()
print("mvec2d_e.new(2,7)*mvec2d.new(1,2,-5,12):",mvec2d_e.new(2,7)*mvec2d.new(1,2,-5,12))
print("mvec2d_e.new(2,7)*vec2d.new(2,-5):",mvec2d_e.new(2,7)*vec2d.new(2,-5))
print("mvec2d.new(1,2,-5,12)*mvec2d_e.new(2,7):",mvec2d.new(1,2,-5,12)*mvec2d_e.new(2,7))
print("vec2d.new(2,-5)*mvec2d_e.new(2,7):",vec2d.new(2,-5)*mvec2d_e.new(2,7))
print()
print("mvec2d_e.new(1,12)*mvec2d_e.new(2,7):",mvec2d_e.new(1,12)*mvec2d_e.new(2,7))
print()
 -- 2d rotations

for i = 0,10
do
    local theta = math.pi/10.0*i
    print("theta:",theta, "exp(I_2d,theta):",exp(I_2d,theta))
end
print()

for i = 0,10
do
    local theta = math.pi/10.0*i
    print("theta:",theta, "rotor(I_2d,theta):",rotor(I_2d,theta),"rotated:", rotate(vec2d.new(1,0),rotor(I_2d,theta)))
end
print()