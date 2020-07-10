function driverRedMaxBDF1(sceneID,batch)
% driverRedMaxBDF1 Reference implementation of the RedMax algorithm
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
	driverRedMaxBDF1(sceneID,true)
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
scene.plotEnergies(1);

end

%%
function simLoop(scene)

jroot = scene.joints{1};
h = scene.h;
nsteps = ceil(scene.tspan(2)/h);

% Integrate
for k = 0 : nsteps
	% Save old state
	[q0,qdot0] = jroot.getQ();
	jroot.setQ0(q0,qdot0);
	
	% Compute new state
	q1 = q0 + h*qdot0; % initial guess
	q1 = newton(@(q1)evalBDF1(q1,scene),q1);
	qdot1 = (q1 - q0)/h;
	
	% Save new state
	jroot.setQ(q1,qdot1);
	
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
function [g,G] = evalBDF1(q1,scene)
h = scene.h;
nr = redmax.Scene.countR();
jroot = scene.joints{1};

% Value from last time step
[q0,qdot0] = jroot.getQ0();

% New values
qdot1 = (q1 - q0)/h;
jroot.setQ(q1,qdot1);
jroot.update();
[Mr,dMrdq,f,K,D] = computeValues(scene);

% BDF1: Newton solve for g(q)=0 with Jacobian G = dg/dq
h2 = h*h;
dqtmp = q1 - q0 - h*qdot0;
g = Mr*dqtmp - h2*f;
G = Mr - h*D - h2*K;
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
