%%
function testRigid

h = 5e-3;
tEnd = 2;
drawHz = 60;
grav = [0 0 -1]' * 98;
rhoV = 1e0; % volume density for rigid bodies

rigids = [];
rigids(end+1).whd = [10 1 1]';
rigids(end).E = eye(4);
rigids(end).phi = [0 3 0 0 0 90]';

% Rigid body mass
nrbs = length(rigids);
for i = 1 : nrbs
	% The rigid frame is attached to the center of mass of the body, with
	% the axes aligned to the inertial frame. Also, we use velocities in
	% body space, so that the mass matrix becomes diagonal.
	whd = rigids(i).whd;
	mass = rhoV * prod(whd);
	rigids(i).m(1) = (1/12) * mass * whd([2,3])' * whd([2,3]);
	rigids(i).m(2) = (1/12) * mass * whd([3,1])' * whd([3,1]);
	rigids(i).m(3) = (1/12) * mass * whd([1,2])' * whd([1,2]);
	rigids(i).m(4) = mass;
	rigids(i).m(5) = mass;
	rigids(i).m(6) = mass;
	fprintf('rigid%d: %.2f g\n', i, mass);
end

% System matrix indices
n = 0;
for i = 1 : nrbs
	% Since we're using maximal coordinates, each rigid body gets 6 degrees
	% of freedom
	rigids(i).i = n + (1:6); %#ok<*AGROW>
	n = n + 6;
end

% Run
t0 = -inf;
for t = 0 : h : tEnd
	% Draw scene
	if t - t0 > 1 / drawHz
		drawTest(t,rigids);
		t0 = t;
	end
	
	% Fill M, f, v
	M = zeros(n);
	f = zeros(n,1);
	v = zeros(n,1);
	for i = 1 : nrbs
		rb = rigids(i);
		R = rb.E(1:3,1:3);
		ir = rb.i;
		phi = rb.phi;
		v(ir) = phi;
		Mr = diag(rb.m);
		M(ir,ir) = Mr; % mass matrix
		fcor = se3.ad(phi)' * Mr * phi; % coriolis force
		fext = [zeros(3,1); rb.m(4) * R' * grav]; % gravity
		f(ir) = fcor + fext;
	end

	% Update velocities and positions
	v = M\(M*v + h*f);
	for i = 1 : nrbs
		rigids(i).phi = v(rigids(i).i);
		rigids(i).E = rigids(i).E*se3.exp(h*rigids(i).phi);
	end

end

end

%%
function drawTest(t,rigids)

if t == 0
	clf;
	xlabel('X');
	ylabel('Y');
	zlabel('Z');
	axis equal;
	axis(20*[-1 1 -1 1 -0 2]);
	grid on;
	view([0,0]);
end
cla;
hold on;

% Draw rigid bodies
for i = 1 : length(rigids)
	se3.drawCuboid(rigids(i).E,rigids(i).whd);
end

str = sprintf('t = %.4f', t);
title(str);
drawnow;

end
