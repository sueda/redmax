function testRedMax(sceneID)
% testRedMax Reference implementation of the Red/Max algorithm
% 
% sceneID: What scene to run.
%          0: Simple serial chain
%          1: Different revolute axes
%          2: Branching

if nargin < 1
	sceneID = 0;
end

scene = testRedMaxScenes(sceneID);
scene.init();

% Draw initial scene
if scene.drawHz ~= 0
	scene.draw(0);
end
fprintf('(%d) ''%s'': tspan = [%.1f %.1f]: nr=%d, nm=%d\n',...
	sceneID,scene.name,scene.tspan,...
	redmax.Scene.countR(),redmax.Scene.countM());

% Integrate
tic
[t,y] = euler(scene);
toc
fprintf('%d steps\n',length(t));

% Draw
dt = scene.tspan(2) - scene.tspan(1);
if scene.drawHz ~= 0 && dt > 0
	for tk = scene.tspan(1) : 1/scene.drawHz : scene.tspan(2)
		[k,s] = searchTime(t,tk);
		ys = (1-s)*y(k,:) + s*y(k+1,:);
		scene.joints{1}.scatterDofs(ys);
		scene.draw(tk);
	end
end

end

%% 
function [i,s] = searchTime(t,ti)
% Finds the index of the time interval around ti
tis = find(t < ti);
if isempty(tis)
	% Beginning of time
	i = 1;
	s = 0;
else
	i = tis(end);
	if i == length(t)
		% End of time
		i = i - 1;
		s = 1;
	else
		% Somewhere in the middle
		t0 = t(i);
		t1 = t(i + 1);
		s = (ti - t0)/(t1 - t0);
	end
end
end

%%
function [t,y] = euler(scene)
nr = redmax.Scene.countR();
body0 = scene.bodies{1};
joint0 = scene.joints{1};
% Allocate t and y
% t(k) is the time at step k
% y(k,:) = [q;qdot]' is the position and velocity at step k
h = scene.hEuler;
t = scene.tspan(1) : h : scene.tspan(2);
y = zeros(length(t),2*nr);
y1 = joint0.gatherDofs();
y(1,:) = y1;
% Integrate
for k = 2 : length(t)
	% Mm: maximal mass matrix
	% fm: maximal force vector (eg gravity)
	% Km: maximal stiffness matrix (eg springs, unused in matlab-simple)
	% Dm: maximal damping matrix (eg body damping)
	% fr: reduced force vector (eg external joint torques)
	% Kr: reduced stiffness matrix (eg joint stiffness)
	% Dr: reduced damping matrix (eg joint damping)
	% J, Jdot: Jacobian and its time derivative
	[Mm,fm] = body0.computeMassGrav(scene.grav);
	[fm,Dm] = body0.computeForceDamping(fm);
	[fr,Kr] = joint0.computeForceStiffness();
	[~,Dr] = joint0.computeForceDamping();
	[J,Jdot] = joint0.computeJacobian();
	q0 = y(k-1,1:nr)';
	qdot0 = y(k-1,nr+1:end)';
	Mr = J'*Mm*J;
	Mr = 0.5*(Mr + Mr');
	frtilde = Mr*qdot0 + h*(J'*(fm - Mm*Jdot*qdot0) + fr);
	Mrtilde = Mr + J'*h*Dm*J + h*Dr - h*h*Kr;
	qdot1 = Mrtilde\frtilde;
	qddot = (qdot1 - qdot0)/h;
	q1 = q0 + h*qdot1;
	yk = [q1;qdot1];
	ydotk = [qdot1;qddot];
	joint0.scatterDofs(yk);
	joint0.scatterDDofs(ydotk);
	y(k,:) = yk;
end
end
