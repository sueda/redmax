function testRedMax(itype,sceneID,drawScene,plotH)
% testRedMax Reference implementation of the Red/Max algorithm
% 
% itype: integrator type
%        1 = recursive O(n) dynamics with ode45
%        2 = Red/Max with ode45
%        3 = Red/Max with Euler
% sceneID: What scene to run. See below for the list.
% drawScene: Whether to draw the scene or not.
% plotH: Whether to plot the total energy over time.

global RECURS_ODE45 REDMAX_ODE45 REDMAX_EULER

RECURS_ODE45 = 1;
REDMAX_ODE45 = 2;
REDMAX_EULER = 3;

if nargin < 4
	plotH = false;
end
if nargin < 3
	drawScene = true;
end
if nargin < 2
	sceneID = 1;
end
if nargin < 1
	itype = 1;
end

% To run in batch mode, run the following:
% clear; clc;
% for itype = 1 : 3
% 	for sceneID = 0 : 34
% 		testRedMax(itype,sceneID,false,false);
% 	end
% end

% Scenes
%  0: Simple serial chain
%  1: Different revolute axes
%  2: Branching
%  3: Spherical joint
%  4: Loop
%  5: Joint torque
%  6: Joint limit
%  7: Equality constrained angles
%  8: Equality constrained angles with a loop
%  9: Hybrid dynamics
% 10: External world force
% 11: Joint stiffness and damping
% 12: Mass-springs
% 13: Maximal hybrid dynamics
% 14: Universal joint
% 15: Prismatic joint
% 16: Planar joint
% 17: Translational joint
% 18: Free joint
% 19: Composite joint
% 20: Reduced/maximal hybrid dynamics
% 21: Spline curve joint
% 22: Spline surface joint
% 23: Point-to-point spring
% 24: Spring damper
% 25: Composite body
% 26: OBJ body
% 27: Internal friction revolute
% 28: Internal friction spherical
% 29: Internal friction prismatic
% 30: Internal friction 4-bar linkage
% 31: External friction
% 32: Prescribed joint via maximal constraint
% 33: External friction 4-bar linkage
% 34: Gears
scene = testRedMaxScenes(itype,sceneID);

% Initialize
scene.init();

% Draw initial scene
if ~drawScene
	scene.drawHz = 0;
end
if scene.drawHz ~= 0
	scene.draw(0);
end
if itype == RECURS_ODE45
	istr = 'RECURS_ODE45';
elseif itype == REDMAX_ODE45
	istr = 'REDMAX_ODE45';
elseif itype == REDMAX_EULER
	istr = 'REDMAX_EULER';
end
fprintf('(%d) ''%s'' with %s: tspan = [%.1f %.1f]: nr=%d, nm=%d\n',...
	sceneID,scene.name,istr,scene.tspan,...
	redmax.Scene.countR(),redmax.Scene.countM());

% Solver options
global qpOpts fminconOpts;
odeOpts = odeset(...
	'RelTol',1e-4,...
	'AbsTol',1e-6);
fminconOpts = optimoptions(@fmincon,...
	'Display','off');
qpOpts = optimoptions(@quadprog,...
	'Display','off');

% Integrate
y0 = scene.joints{1}.gatherDofs();
y0 = scene.deformables{1}.gatherDofs(y0);
tic
if diff(scene.tspan) > 0
	if itype == REDMAX_EULER
		if scene.fric
			[t,y] = eulerFric(scene);
		else
			[t,y] = euler(scene);
		end
	else
		if itype == RECURS_ODE45
			[t,y] = ode45(@(t,y)recursFcn(t,y,scene),scene.tspan,y0,odeOpts);
		elseif itype == REDMAX_ODE45
			[t,y] = ode45(@(t,y)redmaxFcn(t,y,scene),scene.tspan,y0,odeOpts);
		end
	end
else
	t = 0;
	y = y0';
end
toc
fprintf('%d steps\n',length(t));

% Draw
dt = scene.tspan(2) - scene.tspan(1);
if scene.drawHz ~= 0 && dt > 0
	for tk = scene.tspan(1) : 1/scene.drawHz : scene.tspan(2)
		[k,s] = searchTime(t,tk);
		ys = (1-s)*y(k,:) + s*y(k+1,:);
		scene.joints{1}.scatterDofs(ys);
		scene.deformables{1}.scatterDofs(ys);
		scene.draw(tk);
	end
end

% Compute energies
if plotH
	T = zeros(size(t));
	V = zeros(size(t));
	for k = 1 : length(t)
		yk = y(k,:);
		scene.joints{1}.scatterDofs(yk);
		scene.deformables{1}.scatterDofs(yk);
		[T(k),V(k)] = scene.joints{1}.computeEnergies(scene.grav);
		[T(k),V(k)] = scene.deformables{1}.computeEnergies(scene.grav,T(k),V(k));
		V(k) = scene.springs{1}.computeEnergy(V(k));
	end
else
	T = zeros(1,2);
	V = zeros(1,2);
	for k = [1,length(t)]
		yk = y(k,:);
		scene.joints{1}.scatterDofs(yk);
		scene.deformables{1}.scatterDofs(yk);
		[T(k),V(k)] = scene.joints{1}.computeEnergies(scene.grav);
		[T(k),V(k)] = scene.deformables{1}.computeEnergies(scene.grav,T(k),V(k));
		V(k) = scene.springs{1}.computeEnergy(V(k));
	end
end
V = V - V(1);
H = T + V;
dH = H(end) - scene.Hexpected(itype);
if abs(dH) > 1e-2 && dt > 0
	fprintf(2,'### FAIL: %.16f ###\n',H(end));
else
	fprintf(1,'### PASS ###\n');
end
if plotH
	clf;
	plot(t,T,'.-',t,V,'.-',t,T+V,'.-');
	xlabel('time');
	ylabel('energy');
	legend('T','V','T+V');
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
function ydot = recursFcn(t,y,scene)
global qpOpts;
nr = redmax.Scene.countR();
nem = redmax.Scene.countEM();
ner = redmax.Scene.countER();
nim = redmax.Scene.countIM();
nir = redmax.Scene.countIR();
ne = nem + ner;
ni = nim + nir;
body0 = scene.bodies{1};
joint0 = scene.joints{1};
joint1 = scene.joints{end};
spring0 = scene.springs{1};
constraint0 = scene.constraints{1};
joint0.scatterDofs(y);
baum = scene.baumgarte;
scene.sceneFcn(t,scene);
body0.computeMassGrav(scene.grav); % clears external spring
body0.computeForceDamping(); % adds to external spring
spring0.computeForceStiffnessDamping(); % adds to external spring
joint0.rhdPass1();
joint1.rhdPass2(scene.grav);
joint0.rhdPass3();
ydot = joint0.gatherDDofs();
% Handling constraints.
% The constraint matrix A contains inequalities C and equalities G:
%     [Cm*J]
% A = [ Cr ]
%     [Gm*J]
%     [ Gr ]
%        [Cmdot*J + Cm*Jdot]
% Adot = [      Crdot      ]
%        [Gmdot*J + Gm*Jdot]
%        [      Grdot      ]
% We're going to compute A*inv(M)*A' in linear time. M here is the reduced
% mass matrix implicitly computed by the recursive formulation, and 
% W = inv(M). So A*W*A' is actually
%    [Cm*J]
%    [ Cr ]*W*[J'*Cm' Cr' J'*Gm' Gr']
%    [Gm*J]
%    [ Gr ]
% To ensure linear time, we cannot form J or W. Instead, we must multiply
% by J or J' and solve by M, which are all linear time.
if ni > 0
	% Check for active inequality constraints
	[Cm,Cmdot,cm] = constraint0.computeJacIneqM();
	[Cr,Crdot,cr] = constraint0.computeJacIneqR();
	[rowsM,rowsR] = constraint0.getActiveList();
	nim = length(rowsM);
	nir = length(rowsR);
	ni = nim + nir;
	if ni > 0
		Cm = Cm(rowsM,:);
		Cmdot = Cmdot(rowsM,:);
		Cr = Cr(rowsR,:);
		Crdot = Crdot(rowsR,:);
	end
end
if ne > 0
	[Gm,Gmdot,gm] = constraint0.computeJacEqM();
	[Gr,Grdot,gr] = constraint0.computeJacEqR();
end
na = ne + ni; % Total number of constraints
if na > 0
	joint1.computeMinvProdInit();
	AWA = zeros(na);
	ns = [nim,nim+nir,nim+nir+nem,nim+nir+nem+ner];
	for i = 1 :na
		if i <= ns(1)
			% Ai is a row of Cm*J or equivalently a col of J'*Cm'
			Ai = joint1.computeJacTransProd(Cm(i,:)');
		elseif i <= ns(2)
			% Ai is a row of Cr
			Ai = Cr(i-ns(1),:)';
		elseif i <= ns(3)
			% Ai is a row of Gm*J or equivalently a col of J'*Gm'
			Ai = joint1.computeJacTransProd(Gm(i-ns(2),:)');
		else
			% Ai is a row of Gr
			Ai = Gr(i-ns(3),:)';
		end
		joint1.computeMinvProd2(Ai);
		WAi = joint0.computeMinvProd3(Ai);
		if ni > 0
			AWA(      1:ns(1),i) = Cm*joint0.computeJacProd(WAi);
			AWA(ns(1)+1:ns(2),i) = Cr*WAi;
		end
		if ne > 0
			AWA(ns(2)+1:ns(3),i) = Gm*joint0.computeJacProd(WAi);
			AWA(ns(3)+1:ns(4),i) = Gr*WAi;
		end
	end
	qdot = ydot(1:nr);
	qddot = ydot(nr+1:end);
	Jqddot = joint0.computeJacProd(qddot);
	[Jqdot,Jdotqdot] = joint0.computeJacProd(qdot);
	if ni > 0
		% \dot{A*qdot} = [Cm*J]*qddot + [Cmdot*J + Cm*Jdot]*qdot
		%                [ Cr ]         [      Crdot      ]
		dCqdot = [Cm*Jqddot; Cr*qddot] + [Cmdot*Jqdot + Cm*Jdotqdot; Crdot*qdot];
		Cqdot = [Cm*Jqdot; Cr*qdot];
		c = [cm; cr];
	end
	if ne > 0
		% \dot{A*qdot} = [Gm*J]*qddot + [Gmdot*J + Gm*Jdot]*qdot
		%                [ Gr ]         [      Grdot      ]
		dGqdot = [Gm*Jqddot; Gr*qddot] + [Gmdot*Jqdot + Gm*Jdotqdot; Grdot*qdot];
		Gqdot = [Gm*Jqdot; Gr*qdot];
		g = [gm; gr];
	end
	if ne > 0 && ni == 0
		rhs = -dGqdot - 2*baum(1)*Gqdot - baum(2)^2*g;
		l = AWA\rhs;
	elseif ne == 0 && ni > 0
		rhs = dCqdot + 2*baum(1)*Cqdot + baum(2)^2*c;
		l = quadprog(AWA,rhs,[],[],[],[],[],zeros(ni,1),[],qpOpts);
	else
		dAqdot = dCqdot + dGqdot;
		Aqdot = Cqdot + Gqdot;
		a = c + g;
		rhs = dAqdot + 2*baum(1)*Aqdot + baum(2)^2*a;
		l = quadprog(AWA,rhs,[],[],[],[],[],[zeros(ni,1);inf(ne,1)],[],qpOpts);
	end
	% constraint torque
	% taucon = A'*l
	%                                  [lim]
	%        = [J'*Cm' Cr' J'*Gm' Gr']*[lir]
	%                                  [lem]
	%                                  [ler]
	%        = J'*Cm'*lim + Cr'*lir + J'*Gm'*lem + Gr'*ler
	taucon = zeros(nr,1);
	if ni > 0
		taucon = taucon +...
			joint1.computeJacTransProd(Cm'*l(      1:ns(1),1)) +...
			                           Cr'*l(ns(1)+1:ns(2),1);
	end
	if ne > 0
		taucon = taucon + ...
			joint1.computeJacTransProd(Gm'*l(ns(2)+1:ns(3),1)) +...
			                           Gr'*l(ns(3)+1:ns(4),1);
	end
	joint0.scatterTauCon(taucon);
	joint1.rhdPass2(scene.grav);
	joint0.rhdPass3();
	joint0.scatterTauCon(); % reset
	ydot = joint0.gatherDDofs();
end
end

%%
function ydot = redmaxFcn(t,y,scene)
global qpOpts;
nr = redmax.Scene.countR();
nem = redmax.Scene.countEM();
ner = redmax.Scene.countER();
nim = redmax.Scene.countIM();
nir = redmax.Scene.countIR();
ne = nem + ner;
ni = nim + nir;
body0 = scene.bodies{1};
joint0 = scene.joints{1};
deformable0 = scene.deformables{1};
spring0 = scene.springs{1};
constraint0 = scene.constraints{1};
joint0.scatterDofs(y);
deformable0.scatterDofs(y);
baum = scene.baumgarte;
scene.sceneFcn(t,scene);
[Mm,fm] = body0.computeMassGrav(scene.grav);
[fm,Dm] = body0.computeForceDamping(fm);
[Mm,fm] = deformable0.computeMassGrav(scene.grav,Mm,fm);
[fr,Kr] = joint0.computeForceStiffness(); %#ok<ASGLU>
[fr,Dr] = joint0.computeForceDamping(fr); %#ok<ASGLU>
[J,Jdot] = joint0.computeJacobian();
[J,Jdot] = deformable0.computeJacobian(J,Jdot);
[fm,Km,Dm] = spring0.computeForceStiffnessDamping(fm,[],Dm); %#ok<ASGLU>
qdot = y(nr+1:end);
Mr = J'*Mm*J; % Can't apply damping and stiffness matrices with ode45
Mr = 0.5*(Mr + Mr');
fr = J'*(fm - Mm*Jdot*qdot) + fr; % explicit forces
if ne > 0
	[Gm,Gmdot,gm,gmdot,gmddot] = constraint0.computeJacEqM();
	[Gr,Grdot,gr,grdot,grddot] = constraint0.computeJacEqR();
	G = [Gm*J; Gr];
	Gdot = [Gm*Jdot + Gmdot*J; Grdot];
	g = [gm; gr];
	gdot = [gmdot; grdot];
	gddot = [gmddot; grddot];
	rhsG = -(Gdot*qdot + gddot) - 2*baum(1)*(G*qdot + gdot) - baum(2)^2*g;
end
if ni > 0
	% Check for active inequality constraints
	[Cm,Cmdot,cm] = constraint0.computeJacIneqM();
	[Cr,Crdot,cr] = constraint0.computeJacIneqR();
	[rowsM,rowsR] = constraint0.getActiveList();
	nim = length(rowsM);
	nir = length(rowsR);
	ni = nim + nir;
	if ni > 0
		Cm = Cm(rowsM,:);
		Cmdot = Cmdot(rowsM,:);
		Cr = Cr(rowsR,:);
		Crdot = Crdot(rowsR,:);
		C = [Cm*J; Cr];
		Cdot = [Cm*Jdot + Cmdot*J; Crdot];
		c = [cm; cr];
		rhsC = -Cdot*qdot - 2*baum(1)*C*qdot - baum(2)^2*c;
	end
end
if ne == 0 && ni == 0
	% No constraints
	qddot = Mr\fr;
elseif ne > 0 && ni == 0
	% Just equality
	LHS = [Mr G'; G zeros(ne)];
	rhs = [fr; rhsG];
	sol = LHS\rhs;
	qddot = sol(1:nr);
	l = sol(nr+1:end);
	constraint0.scatterForceEqM(Gm',l(1:nem,1));
	constraint0.scatterForceEqR(Gr',l(nem+1:end,1));
elseif ne == 0 && ni > 0
	% Just inequality
	qddot = quadprog(Mr,-fr,C,rhsC,[],[],[],[],[],qpOpts);
else
	% Both equality and inequality
	qddot = quadprog(Mr,-fr,C,rhsC,G,rhsG,[],[],[],qpOpts);
end
ydot = [qdot;qddot];
joint0.scatterDofs(y);
joint0.scatterDDofs(ydot);
deformable0.scatterDofs(y);
deformable0.scatterDDofs(ydot);
end

%%
function [t,y] = euler(scene)
global qpOpts;
nr = redmax.Scene.countR();
nem = redmax.Scene.countEM();
ner = redmax.Scene.countER();
ne = nem + ner;
body0 = scene.bodies{1};
joint0 = scene.joints{1};
deformable0 = scene.deformables{1};
spring0 = scene.springs{1};
constraint0 = scene.constraints{1};
h = scene.hEuler;
t = scene.tspan(1) : h : scene.tspan(2);
y = zeros(length(t),2*nr);
y1 = joint0.gatherDofs();
y1 = deformable0.gatherDofs(y1);
y(1,:) = y1;
baum = scene.baumgarte;
for k = 2 : length(t)
	joint0.reparam();
	yk = joint0.gatherDofs();
	yk = deformable0.gatherDofs(yk);
	y(k-1,:) = yk; % Overwrite with reparameterized state
	scene.sceneFcn(t(k),scene);
	nim = redmax.Scene.countIM();
	nir = redmax.Scene.countIR();
	ni = nim + nir;
	[Mm,fm] = body0.computeMassGrav(scene.grav);
	[tmp,Dm] = body0.computeForceDamping(); % use implicit damping
	[Mm,fm] = deformable0.computeMassGrav(scene.grav,Mm,fm);
	[tmp,Dm] = deformable0.computeForceDamping(tmp,Dm); %#ok<ASGLU>
	[fr,Kr] = joint0.computeForceStiffness();
	[tmp,Dr] = joint0.computeForceDamping(); %#ok<ASGLU>
	[J,Jdot] = joint0.computeJacobian();
	[J,Jdot] = deformable0.computeJacobian(J,Jdot);
	[fm,Km,Dm] = spring0.computeForceStiffnessDamping(fm,[],Dm);
	q0 = y(k-1,1:nr)';
	qdot0 = y(k-1,nr+1:end)';
	Mr = J'*Mm*J;
	Mr = 0.5*(Mr + Mr');
	frtilde = Mr*qdot0 + h*(J'*(fm - Mm*Jdot*qdot0) + fr);
	Mrtilde = Mr + J'*(h*Dm - h*h*Km)*J + h*Dr - h*h*Kr;
	if ne > 0
		[Gm,~,gm,gmdot] = constraint0.computeJacEqM();
		[Gr,~,gr,grdot] = constraint0.computeJacEqR();
		G = [Gm*J; Gr];
		g = [gm; gr];
		gdot = [gmdot; grdot];
		rhsG = -gdot - baum(3)*g;
	end
	if ni > 0
		% Check for active inequality constraints
		[Cm,~,cm,cmdot] = constraint0.computeJacIneqM();
		[Cr,~,cr,crdot] = constraint0.computeJacIneqR();
		[rowsM,rowsR] = constraint0.getActiveList();
		nim = length(rowsM);
		nir = length(rowsR);
		ni = nim + nir;
		if ni > 0
			C = [Cm(rowsM,:)*J; Cr(rowsR,:)];
			c = [cm(rowsM); cr(rowsR)];
			cdot = [cmdot(rowsM); crdot(rowsR)];
			rhsC = -cdot - baum(3)*c;
		end
	end
	if ne == 0 && ni == 0
		% No constraints
		qdot1 = Mrtilde\frtilde;
	elseif ne > 0 && ni == 0
		% Just equality
		LHS = [Mrtilde G'; G zeros(ne)];
		rhs = [frtilde; rhsG];
		sol = LHS\rhs;
		qdot1 = sol(1:nr);
		l = sol(nr+1:end);
		constraint0.scatterForceEqM(Gm',l(1:nem,1)/h);
		constraint0.scatterForceEqR(Gr',l(nem+1:end,1)/h);
	elseif ne == 0 && ni > 0
		% Just inequality
		qdot1 = quadprog(Mrtilde,-frtilde,C,rhsC,[],[],[],[],[],qpOpts);
	else
		% Both equality and inequality
		qdot1 = quadprog(Mrtilde,-frtilde,C,rhsC,G,rhsG,[],[],[],qpOpts);
	end
	qddot = (qdot1 - qdot0)/h;
	q1 = q0 + h*qdot1;
	yk = [q1;qdot1];
	ydotk = [qdot1;qddot];
	joint0.scatterDofs(yk);
	joint0.scatterDDofs(ydotk);
	deformable0.scatterDofs(yk);
	deformable0.scatterDDofs(ydotk);
	y(k,:) = yk;
end
end

%%
function [t,y] = eulerFric(scene)
global qpOpts;
nm = redmax.Scene.countM();
nr = redmax.Scene.countR();
nem = redmax.Scene.countEM();
body0 = scene.bodies{1};
joint0 = scene.joints{1};
joint1 = scene.joints{end};
deformable0 = scene.deformables{1};
spring0 = scene.springs{1};
constraint0 = scene.constraints{1};
h = scene.hEuler;
t = scene.tspan(1) : h : scene.tspan(2);
y = zeros(length(t),2*nr);
y1 = joint0.gatherDofs();
y1 = deformable0.gatherDofs(y1);
y(1,:) = y1;
baum = scene.baumgarte;
fbm = zeros(nm,1); % last friction force
for k = 2 : length(t)
	joint0.reparam();
	%scene.draw(t(k));
	yk = joint0.gatherDofs();
	yk = deformable0.gatherDofs(yk);
	y(k-1,:) = yk; % Overwrite with reparameterized state
	scene.sceneFcn(t(k),scene);
	[Mm,fm] = body0.computeMassGrav(scene.grav);
	[tmp,Dm] = body0.computeForceDamping(); % use implicit damping
	[Mm,fm] = deformable0.computeMassGrav(scene.grav,Mm,fm);
	[tmp,Dm] = deformable0.computeForceDamping(tmp,Dm); %#ok<ASGLU>
	[fr,Kr] = joint0.computeForceStiffness();
	[tmp,Dr] = joint0.computeForceDamping(); %#ok<ASGLU>
	[J,Jdot] = joint0.computeJacobian();
	[J,Jdot] = deformable0.computeJacobian(J,Jdot);
	[fm,Km,Dm] = spring0.computeForceStiffnessDamping(fm,[],Dm);
	q0 = y(k-1,1:nr)';
	qdot0 = y(k-1,nr+1:end)';
	Mr = J'*Mm*J;
	Mr = 0.5*(Mr + Mr');
	frtilde = Mr*qdot0 + h*(J'*(fm - Mm*Jdot*qdot0) + fr);
	Mrtilde = Mr + J'*(h*Dm - h*h*Km)*J + h*Dr - h*h*Kr;
	fmtilde = Mm*J*qdot0 + h*fm;
	Mmtilde = Mm + h*Dm - h*h*Km;
	% Check for active inequality constraints (only maximal)
	[Cm0,~,cm0,cmdot0] = constraint0.computeJacIneqM();
	rowsIM = constraint0.getActiveList();
	constraint0.generateContactsCollision();
	nim = length(rowsIM);
	if nim > 0
		% Keep only active rows
		Cm = Cm0(rowsIM,:);
		rhsC = -cmdot0(rowsIM) - baum(3)*cm0(rowsIM);
		CmJ = Cm*J;
	end
	if nem > 0
		[Gm,~,gm,gmdot] = constraint0.computeJacEqM();
		rhsG = -gmdot - baum(3)*gm;
		GmJ = Gm*J;
	end
	fam0 = zeros(nm,1); % last contact force
	% For clarity, the cases with and without external constraints are 
	% handled separately, even though much of the computations are shared.
	if nem > 0 || nim > 0
		% With external constraints
		T = joint0.computeTangentMatrix();
		T = constraint0.computeTangentMatrix(T);
		nb = size(T,1);
		b = zeros(nb,1);
		for iter = 1 : scene.SPiterMax
			% Compute joint contact force
			% Maximal solve with external constraints
			rhsM = fmtilde + h*fbm;
			if nem > 0 && nim == 0
				LHSm = [Mmtilde Gm'; Gm zeros(nem)];
				rhsm = [rhsM; rhsG];
				solm = LHSm\rhsm;
				vuncb = solm(1:nm);
			elseif nem == 0 && nim > 0
				% The primal problem is:
				%vuncb = quadprog(Mmtilde,-rhsM,Cm,rhsC,[],[],[],[],[],qpOpts);
				% The dual is:
				lclo = zeros(nim,1);
				lchi = inf(nim,1);
				CMC = Cm*(Mmtilde\Cm');
				CMC = 0.5*(CMC + CMC');
				lc = quadprog(CMC,rhsC-Cm*(Mmtilde\rhsM),[],[],[],[],lclo,lchi,[],qpOpts);
				vuncb = Mmtilde\(rhsM - Cm'*lc);
			else
				% The primal problem is:
				%vuncb = quadprog(Mmtilde,-rhsM,Cm,rhsC,Gm,rhsG,[],[],[],qpOpts);
				% The dual is:
				CGm = [Cm; Gm];
				rhsCG = [rhsC; rhsG];
				lcglo = [zeros(nim,1); -inf(nem,1)];
				lcghi = [inf(nim,1); inf(nem,1)];
				CGMGC = CGm*(Mmtilde\CGm');
				CGMGC = 0.5*(CGMGC + CGMGC');
				lcg = quadprog(CGMGC,rhsCG-CGm*(Mmtilde\rhsM),[],[],[],[],lcglo,lcghi,[],qpOpts);
				vuncb = Mmtilde\(rhsM - CGm'*lcg);
			end
			% Reduced solve with external constraints. The difference
			% between vunc and vcon will give us the joint reaction forces.
			% The contact force in the external constraint will need to be
			% dealt with as well.
			if nem > 0 && nim == 0
				LHSr = [Mrtilde GmJ'; GmJ zeros(nem)];
				rhsr = [frtilde + h*J'*fbm; rhsG];
				solr = LHSr\rhsr;
				vconb = J*solr(1:nr);
				lconb = solr(nr+1:end);
				constraint0.scatterForceEqM(Gm',lconb/h);
			elseif nem == 0 && nim > 0
				[solr,~,~,~,lambda] = quadprog(Mrtilde,-(frtilde+h*J'*fbm),CmJ,rhsC,[],[],[],[],[],qpOpts);
				vconb = J*solr;
				lconb = zeros(redmax.Scene.countIM(),1); % active and inactive rows
				lconb(rowsIM) = lambda.ineqlin;
				constraint0.scatterForceIneqM(Cm0',lconb/h);
			else
				[solr,~,~,~,lambda] = quadprog(Mrtilde,-(frtilde+h*J'*fbm),CmJ,rhsC,GmJ,rhsG,[],[],[],qpOpts);
				vconb = J*solr;
				liconb = zeros(redmax.Scene.countIM(),1); % active and inactive rows
				liconb(rowsIM) = lambda.ineqlin;
				leconb = lambda.eqlin;
				constraint0.scatterForceIneqM(Cm0',liconb/h);
				constraint0.scatterForceEqM(Gm',leconb/h);
			end
			constraint0.computeContactMultiplier(h,scene.SPreg);
			fam = Mmtilde*(vconb - vuncb)/h;
			joint1.scatterContactForce(fam);
			joint0.computeContactMultiplier(h,scene.SPreg);
			% Convergence check
			dfam = fam - fam0;
			dfamNormRel = dfam'*(Mmtilde\dfam)/(fam'*(Mmtilde\fam));
			if dfamNormRel < scene.SPconv || iter == scene.SPiterMax
				break;
			end
			fam0 = fam;
			% Compute friction force
			[bl,bu,idx] = joint0.computeFrictionLimits(scene.mu(1),scene.SPathresh);
			[bl,bu,idx] = constraint0.computeFrictionLimits(scene.mu,scene.SPathresh,bl,bu,idx);
			if isempty(idx)
				b = zeros(nb,1);
				fbm = zeros(nm,1);
			else
				T_ = T(idx,:);
				bl_ = bl(idx);
				bu_ = bu(idx);
				H = T_*(Mmtilde\T_');
				H = 0.5*(H + H');
				f = T_*(Mmtilde\(fmtilde + h*fam));
				b_ = quadprog(H,-f,[],[],[],[],bl_,bu_,[],qpOpts);
				b(idx) = b_;
				fbm = -T_'*b_/h;
			end
		end
		if nem > 0 && nim == 0
			LHS = [Mrtilde GmJ'; GmJ zeros(nem)];
			rhs = [frtilde + h*J'*(fam + fbm); rhsG];
			sol = LHS\rhs;
			qdot1 = sol(1:nr);
		elseif nem == 0 && nim > 0
			qdot1 = quadprog(Mrtilde,-(frtilde+h*J'*(fam+fbm)),CmJ,rhsC,[],[],[],[],[],qpOpts);
		else
			qdot1 = quadprog(Mrtilde,-(frtilde+h*J'*(fam+fbm)),CmJ,rhsC,GmJ,rhsG,[],[],[],qpOpts);
		end
	else
		% Without any external constraints
		T = joint0.computeTangentMatrix();
		nb = size(T,1);
		b = zeros(nb,1);
		vunc = Mmtilde\fmtilde;
		vcon = J*(Mrtilde\frtilde);
		for iter = 1 : scene.SPiterMax
			% Compute joint contact force
			vuncb = vunc + h*(Mmtilde\fbm);
			vconb = vcon + h*J*(Mrtilde\(J'*fbm));
			fam = Mmtilde*(vconb - vuncb)/h;
			joint1.scatterContactForce(fam);
			joint0.computeContactMultiplier(h,scene.SPreg);
			% Convergence check
			dfam = fam - fam0;
			dfamNormRel = dfam'*(Mmtilde\dfam)/(fam'*(Mmtilde\fam));
			if dfamNormRel < scene.SPconv || iter == scene.SPiterMax
				break;
			end
			fam0 = fam;
			% Compute friction force
			[bl,bu,idx] = joint0.computeFrictionLimits(scene.mu(1),scene.SPathresh);
			if isempty(idx)
				b = zeros(nb,1);
				fbm = zeros(nm,1);
			else
				T_ = T(idx,:);
				bl_ = bl(idx);
				bu_ = bu(idx);
				H = T_*(Mmtilde\T_');
				H = 0.5*(H + H');
				f = T_*(Mmtilde\(fmtilde + h*fam));
				b_ = quadprog(H,-f,[],[],[],[],bl_,bu_,[],qpOpts);
				b(idx) = b_;
				fbm = -T_'*b_/h;
			end
		end
		qdot1 = Mrtilde\(frtilde + h*J'*(fam + fbm));
	end
	qddot = (qdot1 - qdot0)/h;
	q1 = q0 + h*qdot1;
	yk = [q1;qdot1];
	ydotk = [qdot1;qddot];
	joint0.scatterDofs(yk);
	joint0.scatterDDofs(ydotk);
	deformable0.scatterDofs(yk);
	deformable0.scatterDDofs(ydotk);
	y(k,:) = yk;
end
end
