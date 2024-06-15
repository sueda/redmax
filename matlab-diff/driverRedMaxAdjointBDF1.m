function driverRedMaxAdjointBDF1()
% driverRedMaxAdjointBDF1 Reference implementation of the RedMax algorithm

sceneID = 100;
scene = scenesRedMax(sceneID);
scene.drawHz = 0;
scene.init();
scene.draw();

fprintf('(%d) ''%s'': tEnd=%.1f, nsteps=%d, nr=%d, nm=%d\n',...
	sceneID,scene.name,scene.tEnd,scene.nsteps,...
	redmax.Scene.countR(),redmax.Scene.countM());


% Optimize
pInit = scene.task.p;
opts = optimoptions(@fminunc,...
	'Display','iter-detailed',...
	'SpecifyObjectiveGradient',true,...
	'CheckGradients',false...
	);
tic
p = fminunc(@(p)taskObjective(p,scene),pInit,opts);
toc
fprintf('p = [\n');
disp(p);
fprintf('\b];\n');

% Show the result
scene.drawHz = 15;
scene.reset();
scene.task.p = p;
scene.task.init();
simLoop(scene);

end

%%
function [P,dPdp] = taskObjective(p,scene)
scene.reset();
scene.task.p = p;
scene.task.init();
simLoop(scene);
[P,dPdp] = scene.task.calcFinal();

% Finite difference test
testGrad = false;
if testGrad
	dPdp_ = zeros(1,length(p)); %#ok<UNRCH>
	for i = 1 : length(p)
		p_ = p;
		p_(i) = p_(i) + sqrt(eps);
		scene.reset();
		scene.task.p = p_;
		scene.task.init();
		simLoop(scene);
		P_ = scene.task.calcFinal();
		dPdp_(:,i) = (P_ - P)/sqrt(eps);
	end
	redmax.Scene.printError('dPdp',dPdp_,dPdp);
end
end

%%
function simLoop(scene)

jroot = scene.joints{1};
h = scene.h;
nsteps = scene.nsteps;

% Integrate
for k = 0 : nsteps-1
	% Apply parameters
	scene.task.applyStep();
	
	% Save old state
	[q0,qdot0] = jroot.getQ();
	jroot.setQ0(q0,qdot0);
	
	% Compute new state
	q1 = q0 + h*qdot0; % initial guess
	[q1,GL,GU,Gp,M,f,K,D,J,dJdq] = newton(@(q1)evalBDF1(q1,scene),q1);
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
	scene.saveHistory(GL,GU,Gp,M,f,K,D,J,dJdq);
	scene.draw();
end
%fprintf('%d steps\n',nsteps);

end

%%
function [x,Hl,Hu,Hp,M,f,K,D,J,dJdq] = newton(evalFcn,xInit)
tol = 1e-9;
dxMax = 1e3;
iterMax = 5*length(xInit);
testGrad = false;
x = xInit;
iter = 1;
while true
	[g,H,M,f,K,D,J,dJdq] = evalFcn(x);
	if testGrad
		% Finite difference test
		sqrteps = sqrt(eps); %#ok<UNRCH>
		H_ = zeros(size(H));
		for i = 1 : length(x)
			q_ = x;
			q_(i) = q_(i) + sqrteps;
			g_ = evalFcn(q_);
			H_(:,i) = (g_ - g)/sqrteps;
		end
		redmax.Scene.printError('H',H_,H);
	end
	% dx = -G\g;
	[Hl,Hu,Hp] = lu(H,'vector');
	dx = -(Hu\(Hl\g(Hp)));	
	if norm(dx) > dxMax
		fprintf('Newton diverged\n');
		break;
	end
	% TODO: line search
	x = x + dx;
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
function [g,H,M,f,K,D,J,dJdq] = evalBDF1(q1,scene)
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
	[M,f,dMdq,K,D,J,dJdq] = computeValues(scene);
	g = M*dqtmp - h2*f;
	H = M - h*D - h2*K;
	for i = 1 : nr
		H(:,i) = H(:,i) + dMdq(:,:,i)*dqtmp;
	end
end

end

%%
function [M,f,dMdq,K,D,J,dJdq] = computeValues(scene)
nr = redmax.Scene.countR();
broot = scene.bodies{1};
jroot = scene.joints{1};
froot = scene.forces{1};

qdot = jroot.getQdot();
nm = redmax.Scene.countM();
dJdq = zeros(nm,nr,nr);
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
