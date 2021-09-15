latlon  = [
-11 -174
-3 -175
-7 -176
-1 -177
-9 -178
-12 -179
-6 180
-5 179
-4 178
-8 177
-10 176
-0 175
-12 174
];

lat = latlon(:,1);
lon = latlon(:,2);

colat = 90 - lat;
phi   = (pi/180) .* colat;
theta = (pi/180) .* lon;
x = cos( theta ) .* sin( phi );
y = sin( theta ) .* sin( phi );
z = cos( phi );
rs = sqrt( sum(x)^2 + sum(y)^2 + sum(z)^2 );
xs = sum(x)/rs ;
ys = sum(y)/rs ;
zs = sum(z)/rs ;

mean_lon = atan2( ys, xs ) * (180/pi)       ### theta
mean_lat = 90 - ( acos( zs ) * (180/pi) )   ### phi

n = size(lat);
angular_standard_deviation = acos( rs / n(1) ) * (180/pi)

#### here  kent 
#### theta phi   lon
#### phi   theta lat
####
thetas = acos( zs );
phis   = atan2( ys, xs );

H = [
cos( thetas )*cos( phis ) cos(thetas)*sin(phis) -sin(thetas) 
-sin(phis)                cos(phis)             0
sin(thetas)               sin(thetas)*cos(phis) cos(thetas)
]

T = [
 sum(x.*x) sum(x.*y) sum(x.*z) 
 sum(x.*y) sum(y.*y) sum(y.*z)
 sum(x.*z) sum(y.*z) sum(z.*z)
]

B = H' * (T./n(1)) * H

psi = 0.5 * atan2( (2 * B(1,2)),  (B(1,2) - B(2,2)) )

P = [
cos(psi) -sin(psi) 0
sin(psi)  cos(psi) 0
0         0        1
]

G = H * P

V = G' * (T./n(1)) * G

Q = V(1,1) - V(2,2)

kappa =  (1/(2 - 2*rs - Q)) +  (1/(2 - 2*rs + Q)) 

beta = 0.5 * ( (1/(2 - 2*rs - Q))  -  (1/(2 - 2*rs + Q)) )

