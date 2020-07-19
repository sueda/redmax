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
%    8: Universal joint
%    9: Free3D joint
%    10: Loop
%    11: Free2D with ground (doesn't work with BDF1)
%
%{
% To run in batch mode, run the following:
clear; clc;
for sceneID = 0 : 11
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

fprintf('(%d) ''%s'': tEnd=%.1f, nsteps=%d, nr=%d, nm=%d\n',...
	sceneID,scene.name,scene.tEnd,scene.nsteps,...
	redmax.Scene.countR(),redmax.Scene.countM());

simLoop(scene);
scene.plotEnergies(1);

end

%%
function simLoop(scene)

jroot = scene.joints{1};
h = scene.h;
nsteps = scene.nsteps;

% Integrate
for k = 0 : nsteps-1
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
	scene.k = k + 1;
	
	% End of step
	scene.saveHistory();
	scene.draw();
end
%fprintf('%d steps\n',nsteps);

end

%%
function q = newton(evalFcn,qInit)
tol = 1e-9;
dqMax = 1e3;
iterMax = 5*length(qInit);
iterLsMax = 10;
testGrad = false;
q = qInit;
iter = 1;
while true
	[g,G] = evalFcn(q);
	if testGrad
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
	% Newton direction
	dq = -G\g;
	if norm(dq) > dqMax
		fprintf('Newton diverged\n');
		break;
	end
	% Line search
	alpha = 1;
	g0 = g;
	q0 = q;
	iterLs = 1;
	while true
		q = q0 + alpha*dq;
		g = evalFcn(q);
		if norm(g) < norm(g0)
			break;
		end
		if iterLs >= iterLsMax
			%fprintf('Line search did not converge after %d iterations\n',iterLsMax);
			break;
		end
		alpha = 0.5*alpha;
		iterLs = iterLs + 1;
	end
	if iterLs > 1
		%fprintf('%d line search iters: alpha=%f\n',iterLs,alpha);
	end
	% Convergence check
	if norm(g) < tol
		% Converged
		break;
	end
	if iter >= iterMax
		fprintf('Newton did not converge after %d iterations\n',iterMax);
		break;
	end
	iter = iter + 1;
end
%fprintf('%d\n',iter);
end

%%
function [g,G] = evalBDF1(q1,scene)
h = scene.h;
h2 = h*h;
nr = redmax.Scene.countR();
jroot = scene.joints{1};

% Value from last time step
[q0,qdot0] = jroot.getQ0();
dqtmp = q1 - q0 - h*qdot0;

% New values
qdot1 = (q1 - q0)/h;
jroot.setQ(q1,qdot1);
if nargout == 1
	jroot.update(false);
	[M,f] = computeValues(scene);
	g = M*dqtmp - h2*f;
else
	jroot.update();
	[M,f,dMdq,K,D] = computeValues(scene);
	g = M*dqtmp - h2*f;
	G = M - h*D - h2*K;
	for i = 1 : nr
		G(:,i) = G(:,i) + dMdq(:,:,i)*dqtmp;
	end
end

end

%%
function [M,f,dMdq,K,D] = computeValues(scene)
nr = redmax.Scene.countR();
broot = scene.bodies{1};
jroot = scene.joints{1};
froot = scene.forces{1};

qdot = jroot.getQdot();
if nargout == 2
	[J,Jdot] = jroot.computeJacobian();
	[Mm,fm] = broot.computeMassGrav(scene.grav);
	fm = broot.computeForce(fm);
	fr = jroot.computeForce();
	[fr,fm] = froot.computeValues(fr,fm);
else
	[J,Jdot,dJdq,dJdotdq] = jroot.computeJacobian();
	[Mm,fm,Km,Dm] = broot.computeMassGrav(scene.grav);
	[fm,Km,Dm] = broot.computeForce(fm,Km,Dm);
	[fr,Kr,Dr] = jroot.computeForce();
	[fr,fm,Kr,Km,Dr,Dm] = froot.computeValues(fr,fm,Kr,Km,Dr,Dm);
end

% Inertia
M = J'*Mm*J;

% Forces
fqvv = -J'*Mm*Jdot*qdot;
f = fr + J'*fm + fqvv;

if nargout > 2
	% Derivatives
	dMdq = zeros(nr,nr,nr);
	for i = 1 : nr
		tmp = J'*Mm*dJdq(:,:,i);
		dMdq(:,:,i) = tmp' + tmp;
	end
	
	Kqvv = zeros(nr,nr);
	Dqvv = -J'*Mm*Jdot;
	MmJdotqdot = Mm*Jdot*qdot;
	for i = 1 : nr
		dJdqi = dJdq(:,:,i);
		dJdotdqi = dJdotdq(:,:,i);
		Kqvv(:,i) = -dJdqi'*MmJdotqdot - J'*Mm*dJdotdqi*qdot;
		Dqvv(:,i) = Dqvv(:,i) - J'*Mm*dJdqi*qdot;
	end
	
	K = Kr + J'*Km*J + Kqvv;
	D = Dr + J'*Dm*J + Dqvv;
	for i = 1 : nr
		dJdqi = dJdq(:,:,i);
		K(:,i) = K(:,i) + dJdqi'*fm + J'*Dm*dJdqi*qdot;
	end
end
end
