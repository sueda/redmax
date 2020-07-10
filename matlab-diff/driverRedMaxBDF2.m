function driverRedMaxBDF2(sceneID,batch)
% driverRedMaxBDF2 Reference implementation of the RedMax algorithm
% 
% sceneID: What scene to run.
%    0: Simple serial chain
%    1: Different revolute axes
%    2: Branching
%    3: Prismatic joint
%    4: Planar joint
%    5: Translational joint
%    6: Free2D joint
%    7: Spherical joint
%    8: Loop
%    9: Free3D joint
%
%{
% To run in batch mode, run the following:
clear; clc;
for sceneID = 0 : 9
	driverRedMaxBDF2(sceneID,true)
end
%}

if nargin < 1
	sceneID = 0;
end
if nargin < 2
	batch = false;
end

scene = scenesRedMax(sceneID);
scene.init();
if batch
	scene.drawHz = 0;
	scene.computeH = true;
	scene.plotH = false;
else
	scene.test();
	scene.draw();
end

fprintf('(%d) ''%s'': tspan = [%.1f %.1f]: nr=%d, nm=%d\n',...
	sceneID,scene.name,scene.tspan,...
	redmax.Scene.countR(),redmax.Scene.countM());

simLoop(scene);
scene.plotEnergies(2);

end

%%
function simLoop(scene)

jroot = scene.joints{1};
h = scene.h;
nsteps = ceil(scene.tspan(2)/h);

% Integrate
for k = 0 : nsteps
	if k == 0
		% Take an SDIRK step
		% To compute the values at (1), we take two substeps
		
		% Save old state
		[q0,qdot0] = jroot.getQ();
		jroot.setQ0(q0,qdot0);
		
		% SDIRK2a
		% Values from the last time step are stored in joint.q0 and joint.qdot0
		a = (2-sqrt(2))/2;
		qa = q0 + a*h*qdot0; % initial guess
		qa = newton(@(qa)evalSDIRK2a(qa,scene),qa);
		qdota = (qa - q0)/(a*h);
		
		% Store values from the intermediate time step in joint.q1 and joint.qdot1
		jroot.setQ1(qa,qdota);
		
		% SDIRK2b
		q1 = qa + (1-a)*h*qdota; % initial guess
		q1 = newton(@(q1)evalSDIRK2b(q1,scene),q1);
		qdot1 = (q1 - q0 - (1-a)*h*qdota)/(a*h);
		
		% Save new state
		% Also save to q1 since BDF2 requires q0 and q1
		jroot.setQ(q1,qdot1);
		jroot.setQ1(q0,qdot0);
	else
		% Take a BDF2 step
		% To compute values at (k+1), we need values at (k) and (k-1).
		
		% Save old state: Q1->Q0, Q->Q1
		[q0,qdot0] = jroot.getQ1();
		jroot.setQ0(q0,qdot0);
		[q1,qdot1] = jroot.getQ();
		jroot.setQ1(q1,qdot1);
		
		% BDF2
		q2 = q1 + h*qdot1; % initial guess
		q2 = newton(@(q2)evalBDF2(q2,scene),q2);
		qdot2 = (3/(2*h))*(q2 - (4/3)*q1 + (1/3)*q0);
		
		% Save new state
		jroot.setQ(q2,qdot2);
	end
	
	% Reparameterize if necessary
	jroot.reparam();
	
	% Update time and step
	jroot.update();
	scene.t = scene.t + h;
	scene.k = k;
	
	% End of step
	scene.computeEnergies();
	scene.draw();
end
fprintf('%d steps\n',nsteps);

end

%%
function q = newton(evalFcn,qInit)
tol = 1e-6;
dqMax = 1e3;
iterMax = 5*length(qInit);
testG = false;
q = qInit;
iter = 1;
while true
	[g,G] = evalFcn(q);
	if testG
		% Finite difference test
		sqrteps = sqrt(eps); %#ok<UNRCH>
		G_ = zeros(size(G));
		for i = 1 : length(q)
			q_ = q;
			q_(i) = q_(i) + sqrteps;
			g_ = evalFcn(q_);
			G_(:,i) = (g_ - g)/sqrteps;
		end
		redmax.Scene.printError('G',G_,G);
	end
	if norm(g) < tol
		break;
	end
	dq = -G\g;
	if norm(dq) > dqMax
		fprintf('Newton diverged\n');
		break;
	end
	q = q + dq;
	if norm(dq) < tol
		break;
	end
	if iter >= iterMax
		fprintf('Newton did not converge after %d iterations\n',iterMax);
	end
	iter = iter + 1;
end
%fprintf('%d\n',iter);
end

%%
function [g,G] = evalSDIRK2a(qa,scene)
h = scene.h;
nr = redmax.Scene.countR();
jroot = scene.joints{1};

% Value from last time step
[q0,qdot0] = jroot.getQ0();

% New values
a = (2-sqrt(2))/2;
ah = a*h;
ah2 = ah*ah;
qdota = (qa - q0)/ah;
jroot.setQ(qa,qdota);
jroot.update();
[Mr,dMrdq,f,K,D] = computeValues(scene);

% SDIRK2a: Newton solve for g(q)=0 with Jacobian G = dg/dq
dqtmp = qa - q0 - ah*qdot0;
g = Mr*dqtmp - ah2*f;
G = Mr - ah*D - ah2*K;
for i = 1 : nr
	G(:,i) = G(:,i) + dMrdq(:,:,i)*dqtmp;
end
end

%%
function [g,G] = evalSDIRK2b(q1,scene)
h = scene.h;
nr = redmax.Scene.countR();
jroot = scene.joints{1};

% Values from last time step
[q0,qdot0] = jroot.getQ0();
qdota = jroot.getQdot1();

% New values
a = (2-sqrt(2))/2;
ah = a*h;
ah2 = ah*ah;
qdot1 = (q1 - q0 - (1-a)*h*qdota)/ah;
jroot.setQ(q1,qdot1);
jroot.update();
[Mr,dMrdq,f,K,D] = computeValues(scene);

% SDIRK2b: Newton solve for g(q)=0 with Jacobian G = dg/dq
dqtmp = q1 - q0 - (2*a-1)*h*qdot0 - 2*(1-a)*h*qdota;
g = Mr*dqtmp - ah2*f;
G = Mr - ah*D - ah2*K;
for i = 1 : nr
	G(:,i) = G(:,i) + dMrdq(:,:,i)*dqtmp;
end
end

%%
function [g,G] = evalBDF2(q2,scene)
h = scene.h;
nr = redmax.Scene.countR();
jroot = scene.joints{1};

% Value from last time step
[q0,qdot0] = jroot.getQ0();
[q1,qdot1] = jroot.getQ1();

% New values
qdot2 = (3/(2*h))*(q2 - (4/3)*q1 + (1/3)*q0);
jroot.setQ(q2,qdot2);
jroot.update();
[Mr,dMrdq,f,K,D] = computeValues(scene);

% BDF2: Newton solve for g(q)=0 with Jacobian G = dg/dq
h2 = h*h;
dqtmp = q2 - (4/3)*q1 + (1/3)*q0 - (8/9)*h*qdot1 + (2/9)*h*qdot0;
g = Mr*dqtmp - (4/9)*h2*f;
G = Mr - (2/3)*h*D - (4/9)*h2*K;
for i = 1 : nr
	G(:,i) = G(:,i) + dMrdq(:,:,i)*dqtmp;
end
end

%%
function [Mr,dMrdq,f,K,D] = computeValues(scene)
nr = redmax.Scene.countR();
broot = scene.bodies{1};
jroot = scene.joints{1};
froot = scene.forces{1};

qdot = jroot.getQdot();
[J,dJdq,Jdot,dJdotdq] = jroot.computeJacobian();
[Mm,fm,Km,Dm] = broot.computeMassGrav(scene.grav);
[fm,Km,Dm] = broot.computeForce(fm,Km,Dm);
[fr,Kr,Dr] = jroot.computeForce();
[fr,Kr,Dr,fm,Km,Dm] = froot.computeValues(fr,Kr,Dr,fm,Km,Dm);

% Inertia
Mr = J'*Mm*J;
dMrdq = zeros(nr,nr,nr);
for i = 1 : nr
	tmp = J'*Mm*dJdq(:,:,i);
	dMrdq(:,:,i) = tmp' + tmp;
end

% Quadratic velocity vector
fqvv = -J'*Mm*Jdot*qdot;
Kqvv = zeros(nr,nr);
Dqvv = -J'*Mm*Jdot;
MmJdotqdot = Mm*Jdot*qdot;
for i = 1 : nr
	dJdqi = dJdq(:,:,i);
	dJdotdqi = dJdotdq(:,:,i);
	Kqvv(:,i) = -dJdqi'*MmJdotqdot - J'*Mm*dJdotdqi*qdot;
	Dqvv(:,i) = Dqvv(:,i) - J'*Mm*dJdqi*qdot;
end

% Forces
f = fr + J'*fm + fqvv;
K = Kr + J'*Km*J + Kqvv;
D = Dr + J'*Dm*J + Dqvv;
for i = 1 : nr
	dJdqi = dJdq(:,:,i);
	K(:,i) = K(:,i) + dJdqi'*fm + J'*Dm*dJdqi*qdot;
end
end
