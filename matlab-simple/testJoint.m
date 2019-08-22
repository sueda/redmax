function testJoint

h = 1e-2;
tEnd = 1;
drawHz = 20;
grav = [0 0 -1]' * 98;
rhoV = 1e0; % volume density for rigid bodies

scene = 2;
rigids = [];
joints = [];
switch scene
	case 1
		% Rigid
		rigids(end+1).whd = [10 1 1]';
		rigids(end).E = eye(4);
		rigids(end).phi = [0 5 0 0 0 100]';
	case 2
		% Rigid
		rigids(end+1).whd = [10 1 1]';
		rigids(end).E = eye(4);
		rigids(end).phi = [0 0 0 0 0 0]';
		rigids(end+1).whd = [1 1 10]';
		rigids(end).E = [eye(3), [5 0 -5]'; 0 0 0 1];
		rigids(end).phi = [0 0 0 0 0 0]';
		
		% Joint
		joints(end+1).E = eye(4);
		joints(end).rbs = [0 1];
		joints(end).rows = [1 3 4 5 6];
		joints(end+1).E = [eye(3), [5 0 0]'; 0 0 0 1];
		joints(end).rbs = [1 2];
		joints(end).rows = [1 3 4 5 6];
end

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
mg = 0;
for i = 1 : nrbs
	% Since we're using maximal coordinates, each rigid body gets 6 degrees
	% of freedom
	rigids(i).i = n + (1:6); %#ok<*AGROW>
	n = n + 6;
end

% Initialize joints
njoints = length(joints);
for i = 1 : njoints
	ia = joints(i).rbs(1);
	% E_wa is $^w_aE$, which transforms points from local to world
	if ia == 0
		E_wa = eye(4);
	else
		E_wa = rigids(ia).E;
	end
	E_wb = rigids(joints(i).rbs(2)).E;
	E_wj = joints(i).E;
	E_jw = se3.inv(E_wj);
	% Store where the two bodies are wrt the joint
	joints(i).E0_ja = E_jw * E_wa;
	joints(i).E0_jb = E_jw * E_wb;
	joints(i).mg = mg + (1:length(joints(i).rows));
	% The number of rows in the constraint
	mg = mg + length(joints(i).rows);
end

% Run
t0 = -inf;
for t = 0 : h : tEnd
	% Draw scene
	if t - t0 > 1 / drawHz
		drawTest(t,rigids,joints);
		t0 = t;
	end
	M = zeros(n);
	G = zeros(mg,n);
	g = zeros(mg,1);
	f = zeros(n,1);
	v = zeros(n,1);
	% Rigid bodies
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
	% Joint constraints
	for i = 1 : njoints
		joint = joints(i);
		ia = joint.rbs(1);
		ib = joint.rbs(2);
		rb = rigids(ib);
		if ia == 0
			E_wa = eye(4);
		else
			ra = rigids(ia);
			E_wa = ra.E;
		end
		E_wb = rb.E;
		E_ja = joint.E0_ja; % from a to j
		E_aw = se3.inv(E_wa); % from w to a
		E_jw = E_ja * E_aw; % from w to j
		E_jb = E_jw * E_wb; % from b to j
		Ad_ja = se3.Ad(E_ja);
		Ad_jb = se3.Ad(E_jb);
		rows = joint.rows;
		if ia ~= 0
			G(joint.mg,ra.i) = Ad_ja(rows,:);
		end
		G(joint.mg,rb.i) = -Ad_jb(rows,:);
		% Stabilize just the position
		Ea = E_wa*se3.inv(joint.E0_ja); % joint frame using A
		Eb = E_wb*se3.inv(joint.E0_jb); % joint frame using B
		% The position error only: subtract positions and transform to J
		g(joint.mg(end-2:end)) = E_jw(1:3,1:3)*(Ea(1:3,4) - Eb(1:3,4));
	end
	% Solve KKT system
	ftilde = M*v + h*f;
	LHS = [M G'; G zeros(mg)];
	rhs = [ftilde; -(1/h)*g];
	sol = LHS\rhs;
	v = sol(1:n);
	% Update velocities and positions
	for i = 1 : nrbs
		rigids(i).phi = v(rigids(i).i);
		rigids(i).E = rigids(i).E*se3.exp(h*rigids(i).phi);
	end
end

return

function drawTest(t,rigids,joints)

if t == 0
	clf;
	xlabel('X');
	ylabel('Y');
	zlabel('Z');
	axis equal;
	%axis(20*[-1 1 -1 1 -0 2]);
	grid on;
	view([0,0]);
end
cla;
hold on;

% Draw rigid bodies
for i = 1 : length(rigids)
	se3.drawCuboid(rigids(i).E,rigids(i).whd);
end

% Draw joint
for i = 1 : length(joints)
	joint = joints(i);
	ia = joint.rbs(1);
	ib = joint.rbs(2);
	rb = rigids(ib);
	if ia == 0
		E_wa = eye(4);
	else
		ra = rigids(ia);
		E_wa = ra.E;
	end
	E_wb = rb.E;
	E_aj = se3.inv(joint.E0_ja);
	E_bj = se3.inv(joint.E0_jb);
	se3.drawAxis(E_wa*E_aj,1.0);
	se3.drawAxis(E_wb*E_bj,1.0);
end

str = sprintf('t = %.4f', t);
title(str);
drawnow;

return
