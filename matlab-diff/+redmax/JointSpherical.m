classdef JointSpherical < redmax.Joint
	
	%%
	properties (Constant)
		CHART_XYX = 1;
		CHART_XZX = 2;
		CHART_YZY = 3;
		CHART_YXY = 4;
		CHART_ZXZ = 5;
		CHART_ZYZ = 6;
		CHART_XYZ = 7;
		CHART_XZY = 8;
		CHART_YZX = 9;
		CHART_YXZ = 10;
		CHART_ZXY = 11;
		CHART_ZYX = 12;
	end
	
	%%
	properties
		radius % Radius for display
		chart  % Current coordinate chart
		chart0 % Last chart
		chart1 % Last chart (for BDF2, chart1 is k, and chart0 is k-1)
	end
	
	%%
	methods
		%%
		function this = JointSpherical(parent,body)
			this = this@redmax.Joint(parent,body,3);
			this.radius = 1.0;
			this.chart = redmax.JointSpherical.CHART_XYZ;
		end
		
		%%
		function setGeometry(this,radius)
			this.radius = radius;
		end
		
		%%
		function testEuler(this)
			R = redmax.JointSpherical.getEuler(this.chart,this.q,[0 0 0]');
			q_ = redmax.JointSpherical.getEulerInv(this.chart,R);
			R_ = redmax.JointSpherical.getEuler(this.chart,q_,[0 0 0]');
			redmax.Scene.printError('R',R_,R);
		end
	end
	
	%%
	methods (Access = protected)
		%%
		function setAux0_(this)
			this.chart0 = this.chart;
		end
		
		%%
		function setAux1_(this)
			this.chart1 = this.chart;
		end
		
		%%
		function reparam_(this)
			[R,~,~,~,Told,detTold] = redmax.JointSpherical.getEuler(this.chart,this.q,this.qdot);
			%fprintf('*%s\t%f\t%f %f %f\t%f %f %f\n',this.name,detTold,this.q,this.qdot);
			if abs(detTold) > 0.5
				return;
			end
			
			% Calculate all 12 determinants
			z = zeros(3,1); % qdot not always needed
			% Reconstruct the R matrix from the last time step
			R1 = redmax.JointSpherical.getEuler(this.chart1,this.q1,z);
			detTs = zeros(2,12);
			for k = 1 : 12
				qk = redmax.JointSpherical.getEulerInv(k,R);
				[~,~,~,~,~,detTs(1,k)] = redmax.JointSpherical.getEuler(k,qk,z);
				q1k = redmax.JointSpherical.getEulerInv(k,R1);
				[~,~,~,~,~,detTs(2,k)] = redmax.JointSpherical.getEuler(k,q1k,z);
			end
			% Choose the chart that has the biggest determinant for R and R1.
			detTs(isnan(detTs)) = 0;
			[~,this.chart] = max(min(abs(detTs)));
			fprintf('%s->%s\n',redmax.JointSpherical.getChartName(this.chart0),redmax.JointSpherical.getChartName(this.chart));
			
			% Reparam q
			this.q = redmax.JointSpherical.getEulerInv(this.chart,R);
			% Recompute T
			[~,~,~,~,Tnew] = redmax.JointSpherical.getEuler(this.chart,this.q,z);
			% Reparam qdot
			this.qdot = Tnew\(Told*this.qdot);
			
			% Reparameterize q and qdot from the last step (not needed for
			% BDF1)
			% Assuming BDF2, q is the value just computed, q1 is the
			% previous step, and q0 is from 2 steps before.
			[~,~,~,~,Told] = redmax.JointSpherical.getEuler(this.chart1,this.q1,z);
			this.chart1 = this.chart;
			this.q1 = redmax.JointSpherical.getEulerInv(this.chart1,R1);
			[~,~,~,~,Tnew] = redmax.JointSpherical.getEuler(this.chart1,this.q1,z);
			this.qdot1 = Tnew\(Told*this.qdot1);
		end
		
		%%
		function update_(this)
			[R,dRdq,Rdot,dRdotdq,T,detT,dTdq,Tdot,dTdotdq] = redmax.JointSpherical.getEuler(this.chart,this.q,this.qdot); %#ok<ASGLU>
			%fprintf(' %s\t%f\t%f %f %f\t%f %f %f\n',this.name,detT,this.q,this.qdot);
			n = 3;
			this.Q = eye(4);
			this.Q(1:3,1:3) = R;
			this.A = zeros(6);
			this.A(1:3,1:3) = R;
			this.A(4:6,4:6) = R;
			this.dAdq = zeros(6,6,n);
			this.Adot = zeros(6,6);
			this.Adot(1:3,1:3) = Rdot;
			this.Adot(4:6,4:6) = Rdot;
			this.dAdotdq = zeros(6,6,n);
			this.S = zeros(6,n);
			this.S(1:3,1:3) = T;
			this.dSdq = zeros(6,n,n);
			this.Sdot = zeros(6,n);
			this.Sdot(1:3,1:3) = Tdot;
			this.dSdotdq = zeros(6,n,n);
			for k = 1 : n
				dRdqk = dRdq(:,:,k);
				dRdotdqk = dRdotdq(:,:,k);
				this.dAdq(1:3,1:3,k) = dRdqk;
				this.dAdq(4:6,4:6,k) = dRdqk;
				this.dAdotdq(1:3,1:3,k) = dRdotdqk;
				this.dAdotdq(4:6,4:6,k) = dRdotdqk;
				this.dSdq(1:3,:,k) = dTdq(:,:,k);
				this.dSdotdq(1:3,:,k) = dTdotdq(:,:,k);
			end
		end
		
		%%
		function draw_(this)
			E = this.E_wj;
			r = this.radius;
			S = diag([r,r,r,1]);
			n = 8;
			[X,Y,Z] = sphere(n);
			X = reshape(X,1,(n+1)^2);
			Y = reshape(Y,1,(n+1)^2);
			Z = reshape(Z,1,(n+1)^2);
			XYZ = E*S*[X;Y;Z;ones(1,(n+1)^2)];
			X = reshape(XYZ(1,:),n+1,n+1);
			Y = reshape(XYZ(2,:),n+1,n+1);
			Z = reshape(XYZ(3,:),n+1,n+1);
			surf(X,Y,Z,'FaceColor',this.body.color);
		end
	end
	
	methods (Static)
		%%
		function [R,dRdq,Rdot,dRdotdq,T,detT,dTdq,Tdot,dTdotdq] = getEuler(chart,q,qdot)
			switch chart
				case redmax.JointSpherical.CHART_XYX
					[R,dRdq,Rdot,dRdotdq,T,detT,dTdq,Tdot,dTdotdq] = redmax.JointSpherical.XYX(q,qdot);
				case redmax.JointSpherical.CHART_XZX
					[R,dRdq,Rdot,dRdotdq,T,detT,dTdq,Tdot,dTdotdq] = redmax.JointSpherical.XZX(q,qdot);
				case redmax.JointSpherical.CHART_YZY
					[R,dRdq,Rdot,dRdotdq,T,detT,dTdq,Tdot,dTdotdq] = redmax.JointSpherical.YZY(q,qdot);
				case redmax.JointSpherical.CHART_YXY
					[R,dRdq,Rdot,dRdotdq,T,detT,dTdq,Tdot,dTdotdq] = redmax.JointSpherical.YXY(q,qdot);
				case redmax.JointSpherical.CHART_ZXZ
					[R,dRdq,Rdot,dRdotdq,T,detT,dTdq,Tdot,dTdotdq] = redmax.JointSpherical.ZXZ(q,qdot);
				case redmax.JointSpherical.CHART_ZYZ
					[R,dRdq,Rdot,dRdotdq,T,detT,dTdq,Tdot,dTdotdq] = redmax.JointSpherical.ZYZ(q,qdot);
				case redmax.JointSpherical.CHART_XYZ
					[R,dRdq,Rdot,dRdotdq,T,detT,dTdq,Tdot,dTdotdq] = redmax.JointSpherical.XYZ(q,qdot);
				case redmax.JointSpherical.CHART_XZY
					[R,dRdq,Rdot,dRdotdq,T,detT,dTdq,Tdot,dTdotdq] = redmax.JointSpherical.XZY(q,qdot);
				case redmax.JointSpherical.CHART_YZX
					[R,dRdq,Rdot,dRdotdq,T,detT,dTdq,Tdot,dTdotdq] = redmax.JointSpherical.YZX(q,qdot);
				case redmax.JointSpherical.CHART_YXZ
					[R,dRdq,Rdot,dRdotdq,T,detT,dTdq,Tdot,dTdotdq] = redmax.JointSpherical.YXZ(q,qdot);
				case redmax.JointSpherical.CHART_ZXY
					[R,dRdq,Rdot,dRdotdq,T,detT,dTdq,Tdot,dTdotdq] = redmax.JointSpherical.ZXY(q,qdot);
				case redmax.JointSpherical.CHART_ZYX
					[R,dRdq,Rdot,dRdotdq,T,detT,dTdq,Tdot,dTdotdq] = redmax.JointSpherical.ZYX(q,qdot);
			end
		end
		
		%%
		function q = getEulerInv(chart,R)
			switch chart
				case redmax.JointSpherical.CHART_XYX
					q = redmax.JointSpherical.XYXinv(R);
				case redmax.JointSpherical.CHART_XZX
					q = redmax.JointSpherical.XZXinv(R);
				case redmax.JointSpherical.CHART_YZY
					q = redmax.JointSpherical.YZYinv(R);
				case redmax.JointSpherical.CHART_YXY
					q = redmax.JointSpherical.YXYinv(R);
				case redmax.JointSpherical.CHART_ZXZ
					q = redmax.JointSpherical.ZXZinv(R);
				case redmax.JointSpherical.CHART_ZYZ
					q = redmax.JointSpherical.ZYZinv(R);
				case redmax.JointSpherical.CHART_XYZ
					q = redmax.JointSpherical.XYZinv(R);
				case redmax.JointSpherical.CHART_XZY
					q = redmax.JointSpherical.XZYinv(R);
				case redmax.JointSpherical.CHART_YZX
					q = redmax.JointSpherical.YZXinv(R);
				case redmax.JointSpherical.CHART_YXZ
					q = redmax.JointSpherical.YXZinv(R);
				case redmax.JointSpherical.CHART_ZXY
					q = redmax.JointSpherical.ZXYinv(R);
				case redmax.JointSpherical.CHART_ZYX
					q = redmax.JointSpherical.ZYXinv(R);
			end
		end
		
		%%
		function s = getChartName(chart)
			switch chart
				case redmax.JointSpherical.CHART_XYX
					s = 'XYX';
				case redmax.JointSpherical.CHART_XZX
					s = 'XZX';
				case redmax.JointSpherical.CHART_YZY
					s = 'YZY';
				case redmax.JointSpherical.CHART_YXY
					s = 'YXY';
				case redmax.JointSpherical.CHART_ZXZ
					s = 'ZXZ';
				case redmax.JointSpherical.CHART_ZYZ
					s = 'ZYZ';
				case redmax.JointSpherical.CHART_XYZ
					s = 'XYZ';
				case redmax.JointSpherical.CHART_XZY
					s = 'XZY';
				case redmax.JointSpherical.CHART_YZX
					s = 'YZX';
				case redmax.JointSpherical.CHART_YXZ
					s = 'YXZ';
				case redmax.JointSpherical.CHART_ZXY
					s = 'ZXY';
				case redmax.JointSpherical.CHART_ZYX
					s = 'ZYX';
			end
		end
		
		%%
		function codegen()
			clear;
			q = sym('q',[3,1],'real');
			qdot = sym('qdot',[3,1],'real');
			s = sin(q);
			c = cos(q);
			X1 = [1 0 0; 0 c(1) -s(1); 0 s(1) c(1)];
			X2 = [1 0 0; 0 c(2) -s(2); 0 s(2) c(2)];
			X3 = [1 0 0; 0 c(3) -s(3); 0 s(3) c(3)];
			Y1 = [c(1) 0 s(1); 0 1 0; -s(1) 0 c(1)];
			Y2 = [c(2) 0 s(2); 0 1 0; -s(2) 0 c(2)];
			Y3 = [c(3) 0 s(3); 0 1 0; -s(3) 0 c(3)];
			Z1 = [c(1) -s(1) 0; s(1) c(1) 0; 0 0 1];
			Z2 = [c(2) -s(2) 0; s(2) c(2) 0; 0 0 1];
			Z3 = [c(3) -s(3) 0; s(3) c(3) 0; 0 0 1];
			
			XYX = X1*Y2*X3; XZX = X1*Z2*X3;
			YZY = Y1*Z2*Y3; YXY = Y1*X2*Y3;
			ZXZ = Z1*X2*Z3; ZYZ = Z1*Y2*Z3;
			XYZ = X1*Y2*Z3; XZY = X1*Z2*Y3;
			YZX = Y1*Z2*X3; YXZ = Y1*X2*Z3;
			ZXY = Z1*X2*Y3; ZYX = Z1*Y2*X3;
			
			Rs(1:6)     = { XYX   XZX   YZY   YXY   ZXZ   ZYZ };
			Rs(7:12)    = { XYZ   XZY   YZX   YXZ   ZXY   ZYX };
			names(1:6)  = {'XYX' 'XZX' 'YZY' 'YXY' 'ZXZ' 'ZYZ'};
			names(7:12) = {'XYZ' 'XZY' 'YZX' 'YXZ' 'ZXY' 'ZYX'};
			dRdqs = cell(size(Rs));
			Rdots = cell(size(Rs));
			dRdotdqs = cell(size(Rs));
			Ss = cell(size(Rs));
			detSs = cell(size(Rs));
			dSdqs = cell(size(Rs));
			Sdots = cell(size(Rs));
			dSdotdqs = cell(size(Rs));
			
			for k = 1 : length(Rs)
				R = Rs{k};
				% dRdq
				dRdq = 0*sym('dRdq',[3,3,3],'real');
				for i = 1 : 3
					dRdq(:,:,i) = diff(R,q(i));
				end
				dRdqs{k} = dRdq;
				% Rdot
				Rdot = 0*sym('Rdot',[3,3],'real');
				for i = 1 : 3
					Rdot = Rdot + dRdq(:,:,i)*qdot(i);
				end
				Rdot = simplify(Rdot);
				Rdots{k} = Rdot;
				% dRdotdq
				dRdotdq = 0*sym('dRdotdq',[3,3,3],'real');
				for i = 1 : 3
					dRdotdq(:,:,i) = diff(Rdot,q(i));
				end
				dRdotdqs{k} = dRdotdq;
				% S
				S = 0*sym('S',[3,3],'real');
				for i = 1 : 3
					tmp = simplify(R'*diff(R,q(i)));
					S(1:3,i) = [tmp(3,2); tmp(1,3); tmp(2,1)];
				end
				Ss{k} = S;
				% det(S)
				detS = simplify(det(S));
				detSs{k} = detS;
				% dSdq
				dSdq = 0*sym('dSdq',[3,3,3],'real');
				for i = 1 : 3
					dSdq(:,:,i) = diff(S,q(i));
				end
				dSdqs{k} = dSdq;
				% Sdot
				Sdot = 0*sym('Sdot',[3,3],'real');
				for i = 1 : 3
					Sdot = Sdot + dSdq(:,:,i)*qdot(i);
				end
				Sdots{k} = Sdot;
				% dSdotdq
				dSdotdq = 0*sym('dSdotdq',[3,3,3],'real');
				for i = 1 : 3
					dSdotdq(:,:,i) = diff(Sdot,q(i));
				end
				dSdotdqs{k} = dSdotdq;
				% Output
				filename = sprintf('%s.m',names{k});
				matlabFunction(R,dRdq,Rdot,dRdotdq,...
					S,detS,dSdq,Sdot,dSdotdq,...
					'File',filename,'Optimize',true,'Vars',{q,qdot},...
					'Outputs',{'R','dRdq','Rdot','dRdotdq',...
					'S','detS','dSdq','Sdot','dSdotdq'});
				fprintf('%s\n',names{k});
			end
		end
	end
	
	%%
	methods (Static, Access = private)
		
		%%
		function [R,dRdq,Rdot,dRdotdq,S,detS,dSdq,Sdot,dSdotdq] = XYX(in1,in2)
			q1 = in1(1,:);
			q2 = in1(2,:);
			q3 = in1(3,:);
			qdot1 = in2(1,:);
			qdot2 = in2(2,:);
			qdot3 = in2(3,:);
			t2 = cos(q1);
			t3 = cos(q2);
			t4 = cos(q3);
			t5 = sin(q1);
			t6 = sin(q2);
			t7 = sin(q3);
			t8 = qdot2.*t3;
			t9 = qdot3.*t4;
			t10 = qdot2.*t6;
			t11 = qdot3.*t7;
			t12 = t2.*t4;
			t13 = t3.*t4;
			t14 = t2.*t6;
			t15 = t2.*t7;
			t16 = t4.*t5;
			t17 = t3.*t7;
			t18 = t4.*t6;
			t19 = t5.*t6;
			t20 = t5.*t7;
			t21 = t6.*t7;
			t22 = -t6;
			t23 = -t7;
			t24 = t2.*t8;
			t25 = t4.*t8;
			t26 = t3.*t9;
			t27 = qdot1.*t14;
			t28 = t5.*t8;
			t29 = t7.*t8;
			t30 = t4.*t10;
			t31 = t3.*t11;
			t32 = t6.*t9;
			t33 = qdot1.*t19;
			t36 = -t8;
			t37 = -t9;
			t38 = -t10;
			t39 = -t15;
			t40 = -t16;
			t41 = -t20;
			t42 = -t21;
			t43 = t3.*t12;
			t44 = t3.*t15;
			t45 = t5.*t13;
			t46 = t5.*t17;
			t47 = t10.*t12;
			t48 = t10.*t15;
			t49 = t10.*t16;
			t50 = t10.*t20;
			t51 = t8.*t23;
			t54 = t9.*t22;
			t55 = t10.*t23;
			t56 = t11.*t22;
			t52 = -t30;
			t53 = -t31;
			t57 = -t43;
			t58 = t3.*t39;
			t59 = -t45;
			t60 = -t46;
			t61 = t12.*t38;
			t62 = t15+t45;
			t63 = t16+t44;
			t64 = t27+t28;
			t65 = t29+t32;
			t72 = t41+t43;
			t73 = t25+t56;
			t74 = t26+t55;
			t85 = t51+t54;
			t66 = qdot1.*t62;
			t67 = qdot1.*t63;
			t68 = qdot3.*t62;
			t69 = qdot3.*t63;
			t70 = t12+t60;
			t71 = t20+t57;
			t83 = t39+t59;
			R = reshape([t3,t19,-t14,t21,t70,t63,t18,t83,t72],[3,3]);
			if nargout > 1
				t84 = t40+t58;
				dRdq = reshape([0.0,t14,t19,0.0,t84,t70,0.0,t71,t83,t22,t3.*t5,-t2.*t3,t17,t7.*t19,t14.*t23,t13,t6.*t16,t12.*t22,0.0,0.0,0.0,t18,t83,t72,t42,-t12+t46,t84],[3,3,3]);
			end
			if nargout > 2
				t86 = t52+t53;
				t75 = qdot1.*t70;
				t76 = qdot1.*t71;
				t77 = -t66;
				t78 = -t67;
				t79 = qdot3.*t70;
				t80 = qdot3.*t71;
				t81 = -t68;
				t82 = -t69;
				t87 = -t75;
				t88 = -t79;
				t89 = t50+t78+t81;
				t90 = t61+t77+t82;
				t91 = t49+t76+t88;
				Rdot = reshape([t38,t64,-t24+t33,t65,t89,t75-t80+t15.*t38,t73,t91,t90],[3,3]);
			end
			if nargout > 3
				t92 = t48+t80+t87;
				dRdotdq = reshape([0.0,t24-t33,t64,0.0,t92,t89,0.0,t47+t66+t69,t91,t36,t5.*t38+qdot1.*t2.*t3,t2.*t10+qdot1.*t3.*t5,t74,t8.*t20+t9.*t19+t7.*t27,t7.*t33+t14.*t37+t15.*t36,t86,t8.*t16-t11.*t19+qdot1.*t6.*t12,t11.*t14+t12.*t36+qdot1.*t6.*t16,0.0,0.0,0.0,t73,t91,t90,t85,t67+t68+t20.*t38,t92],[3,3,3]);
			end
			if nargout > 4
				S = reshape([t3,t21,t18,0.0,t4,t23,1.0,0.0,0.0],[3,3]);
			end
			if nargout > 5
				detS = t22;
			end
			if nargout > 6
				dSdq = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t22,t17,t13,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t18,t42,0.0,t23,-t4,0.0,0.0,0.0],[3,3,3]);
			end
			if nargout > 7
				Sdot = reshape([t38,t65,t73,0.0,-t11,t37,0.0,0.0,0.0],[3,3]);
			end
			if nargout > 8
				dSdotdq = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t36,t74,t86,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t73,t85,0.0,t37,t11,0.0,0.0,0.0],[3,3,3]);
			end
		end
		
		%%
		function [R,dRdq,Rdot,dRdotdq,S,detS,dSdq,Sdot,dSdotdq] = XZX(in1,in2)
			q1 = in1(1,:);
			q2 = in1(2,:);
			q3 = in1(3,:);
			qdot1 = in2(1,:);
			qdot2 = in2(2,:);
			qdot3 = in2(3,:);
			t2 = cos(q1);
			t3 = cos(q2);
			t4 = cos(q3);
			t5 = sin(q1);
			t6 = sin(q2);
			t7 = sin(q3);
			t8 = qdot2.*t3;
			t9 = qdot3.*t4;
			t10 = qdot2.*t6;
			t11 = qdot3.*t7;
			t12 = t2.*t4;
			t13 = t3.*t4;
			t14 = t2.*t6;
			t15 = t2.*t7;
			t16 = t4.*t5;
			t17 = t3.*t7;
			t18 = t4.*t6;
			t19 = t5.*t6;
			t20 = t5.*t7;
			t21 = t6.*t7;
			t22 = -t6;
			t23 = t2.*t8;
			t24 = t4.*t8;
			t25 = t3.*t9;
			t26 = qdot1.*t14;
			t27 = t5.*t8;
			t28 = t7.*t8;
			t29 = t4.*t10;
			t30 = t3.*t11;
			t31 = t6.*t9;
			t32 = qdot1.*t19;
			t33 = t7.*t10;
			t34 = t6.*t11;
			t35 = -t8;
			t36 = -t10;
			t37 = -t11;
			t38 = -t13;
			t39 = -t15;
			t40 = -t16;
			t41 = -t18;
			t42 = -t20;
			t43 = t3.*t12;
			t44 = t3.*t15;
			t45 = t5.*t13;
			t46 = t5.*t17;
			t47 = t10.*t12;
			t48 = t10.*t15;
			t49 = t10.*t16;
			t51 = t10.*t20;
			t54 = t11.*t22;
			t50 = -t24;
			t52 = -t32;
			t53 = -t33;
			t55 = -t43;
			t56 = t3.*t39;
			t57 = t5.*t38;
			t58 = -t46;
			t59 = t12.*t36;
			t60 = t15+t45;
			t61 = t16+t44;
			t62 = t28+t31;
			t63 = t29+t30;
			t70 = t42+t43;
			t72 = t24+t54;
			t64 = qdot1.*t60;
			t65 = qdot1.*t61;
			t66 = qdot3.*t60;
			t67 = qdot3.*t61;
			t68 = t12+t58;
			t69 = t20+t55;
			t71 = t23+t52;
			t73 = t25+t53;
			t74 = t34+t50;
			t83 = t39+t57;
			t84 = t40+t56;
			R = reshape([t3,t14,t19,t41,t70,t60,t21,t84,t68],[3,3]);
			if nargout > 1
				dRdq = reshape([0.0,-t19,t14,0.0,t83,t70,0.0,-t12+t46,t84,t22,t2.*t3,t3.*t5,t38,t12.*t22,t16.*t22,t17,t7.*t14,t7.*t19,0.0,0.0,0.0,t21,t84,t68,t18,t69,t83],[3,3,3]);
			end
			if nargout > 2
				t75 = qdot1.*t68;
				t76 = qdot1.*t69;
				t77 = -t64;
				t78 = -t65;
				t79 = qdot3.*t68;
				t80 = qdot3.*t69;
				t81 = -t66;
				t82 = -t67;
				t85 = -t75;
				t86 = -t79;
				t87 = t51+t78+t81;
				t88 = t59+t77+t82;
				t89 = t49+t76+t86;
				t90 = t48+t80+t85;
				Rdot = reshape([t36,t71,t26+t27,t74,t88,-t76+t79+t16.*t36,t62,t90,t87],[3,3]);
			end
			if nargout > 3
				dRdotdq = reshape([0.0,-t26-t27,t71,0.0,t89,t88,0.0,t65+t66+t20.*t36,t90,t35,t2.*t36-qdot1.*t3.*t5,t5.*t36+qdot1.*t2.*t3,t63,t11.*t14+t12.*t35+qdot1.*t6.*t16,t11.*t19+t16.*t35+qdot1.*t12.*t22,t73,t8.*t15+t9.*t14+t7.*t52,t8.*t20+t9.*t19+t7.*t26,0.0,0.0,0.0,t62,t90,t87,t72,t47+t64+t67,t89],[3,3,3]);
			end
			if nargout > 4
				S = reshape([t3,t41,t21,0.0,t7,t4,1.0,0.0,0.0],[3,3]);
			end
			if nargout > 5
				detS = t22;
			end
			if nargout > 6
				dSdq = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t22,t38,t17,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t21,t18,0.0,t4,-t7,0.0,0.0,0.0],[3,3,3]);
			end
			if nargout > 7
				Sdot = reshape([t36,t74,t62,0.0,t9,t37,0.0,0.0,0.0],[3,3]);
			end
			if nargout > 8
				dSdotdq = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t35,t63,t73,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t62,t72,0.0,t37,-t9,0.0,0.0,0.0],[3,3,3]);
			end
		end
		
		%%
		function [R,dRdq,Rdot,dRdotdq,S,detS,dSdq,Sdot,dSdotdq] = YZY(in1,in2)
			q1 = in1(1,:);
			q2 = in1(2,:);
			q3 = in1(3,:);
			qdot1 = in2(1,:);
			qdot2 = in2(2,:);
			qdot3 = in2(3,:);
			t2 = cos(q1);
			t3 = cos(q2);
			t4 = cos(q3);
			t5 = sin(q1);
			t6 = sin(q2);
			t7 = sin(q3);
			t8 = qdot2.*t3;
			t9 = qdot3.*t4;
			t10 = qdot2.*t6;
			t11 = qdot3.*t7;
			t12 = t2.*t4;
			t13 = t3.*t4;
			t14 = t2.*t6;
			t15 = t2.*t7;
			t16 = t4.*t5;
			t17 = t3.*t7;
			t18 = t4.*t6;
			t19 = t5.*t6;
			t20 = t5.*t7;
			t21 = t6.*t7;
			t22 = -t6;
			t23 = -t7;
			t24 = t2.*t8;
			t25 = t4.*t8;
			t26 = t3.*t9;
			t27 = qdot1.*t14;
			t28 = t5.*t8;
			t29 = t7.*t8;
			t30 = t4.*t10;
			t31 = t3.*t11;
			t32 = t6.*t9;
			t33 = qdot1.*t19;
			t36 = -t8;
			t37 = -t9;
			t38 = -t10;
			t39 = -t15;
			t40 = -t16;
			t41 = -t20;
			t42 = -t21;
			t43 = t3.*t12;
			t44 = t3.*t15;
			t45 = t5.*t13;
			t46 = t5.*t17;
			t47 = t10.*t12;
			t48 = t10.*t15;
			t49 = t10.*t16;
			t50 = t10.*t20;
			t51 = t8.*t23;
			t54 = t9.*t22;
			t55 = t10.*t23;
			t56 = t11.*t22;
			t52 = -t30;
			t53 = -t31;
			t57 = -t43;
			t58 = t3.*t39;
			t59 = -t45;
			t60 = -t46;
			t61 = t12.*t38;
			t62 = t15+t45;
			t63 = t16+t44;
			t64 = t27+t28;
			t65 = t29+t32;
			t72 = t41+t43;
			t73 = t25+t56;
			t74 = t26+t55;
			t85 = t51+t54;
			t66 = qdot1.*t62;
			t67 = qdot1.*t63;
			t68 = qdot3.*t62;
			t69 = qdot3.*t63;
			t70 = t12+t60;
			t71 = t20+t57;
			t83 = t39+t59;
			R = reshape([t72,t18,t83,-t14,t3,t19,t63,t21,t70],[3,3]);
			if nargout > 1
				t84 = t40+t58;
				dRdq = reshape([t83,0.0,t71,t19,0.0,t14,t70,0.0,t84,t12.*t22,t13,t6.*t16,-t2.*t3,t22,t3.*t5,t14.*t23,t17,t7.*t19,t84,t42,-t12+t46,0.0,0.0,0.0,t72,t18,t83],[3,3,3]);
			end
			if nargout > 2
				t86 = t52+t53;
				t75 = qdot1.*t70;
				t76 = qdot1.*t71;
				t77 = -t66;
				t78 = -t67;
				t79 = qdot3.*t70;
				t80 = qdot3.*t71;
				t81 = -t68;
				t82 = -t69;
				t87 = -t75;
				t88 = -t79;
				t89 = t50+t78+t81;
				t90 = t61+t77+t82;
				t91 = t49+t76+t88;
				Rdot = reshape([t90,t73,t91,-t24+t33,t38,t64,t75-t80+t15.*t38,t65,t89],[3,3]);
			end
			if nargout > 3
				t92 = t48+t80+t87;
				dRdotdq = reshape([t91,0.0,t47+t66+t69,t64,0.0,t24-t33,t89,0.0,t92,t11.*t14+t12.*t36+qdot1.*t6.*t16,t86,t8.*t16-t11.*t19+qdot1.*t6.*t12,t2.*t10+qdot1.*t3.*t5,t36,t5.*t38+qdot1.*t2.*t3,t7.*t33+t14.*t37+t15.*t36,t74,t8.*t20+t9.*t19+t7.*t27,t92,t85,t67+t68+t20.*t38,0.0,0.0,0.0,t90,t73,t91],[3,3,3]);
			end
			if nargout > 4
				S = reshape([t18,t3,t21,t23,0.0,t4,0.0,1.0,0.0],[3,3]);
			end
			if nargout > 5
				detS = t22;
			end
			if nargout > 6
				dSdq = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t13,t22,t17,0.0,0.0,0.0,0.0,0.0,0.0,t42,0.0,t18,-t4,0.0,t23,0.0,0.0,0.0],[3,3,3]);
			end
			if nargout > 7
				Sdot = reshape([t73,t38,t65,t37,0.0,-t11,0.0,0.0,0.0],[3,3]);
			end
			if nargout > 8
				dSdotdq = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t86,t36,t74,0.0,0.0,0.0,0.0,0.0,0.0,t85,0.0,t73,t11,0.0,t37,0.0,0.0,0.0],[3,3,3]);
			end
		end
		
		%%
		function [R,dRdq,Rdot,dRdotdq,S,detS,dSdq,Sdot,dSdotdq] = YXY(in1,in2)
			q1 = in1(1,:);
			q2 = in1(2,:);
			q3 = in1(3,:);
			qdot1 = in2(1,:);
			qdot2 = in2(2,:);
			qdot3 = in2(3,:);
			t2 = cos(q1);
			t3 = cos(q2);
			t4 = cos(q3);
			t5 = sin(q1);
			t6 = sin(q2);
			t7 = sin(q3);
			t8 = qdot2.*t3;
			t9 = qdot3.*t4;
			t10 = qdot2.*t6;
			t11 = qdot3.*t7;
			t12 = t2.*t4;
			t13 = t3.*t4;
			t14 = t2.*t6;
			t15 = t2.*t7;
			t16 = t4.*t5;
			t17 = t3.*t7;
			t18 = t4.*t6;
			t19 = t5.*t6;
			t20 = t5.*t7;
			t21 = t6.*t7;
			t22 = -t6;
			t23 = t2.*t8;
			t24 = t4.*t8;
			t25 = t3.*t9;
			t26 = qdot1.*t14;
			t27 = t5.*t8;
			t28 = t7.*t8;
			t29 = t4.*t10;
			t30 = t3.*t11;
			t31 = t6.*t9;
			t32 = qdot1.*t19;
			t33 = t7.*t10;
			t34 = t6.*t11;
			t35 = -t8;
			t36 = -t10;
			t37 = -t11;
			t38 = -t13;
			t39 = -t15;
			t40 = -t16;
			t41 = -t18;
			t42 = -t20;
			t43 = t3.*t12;
			t44 = t3.*t15;
			t45 = t5.*t13;
			t46 = t5.*t17;
			t47 = t10.*t12;
			t48 = t10.*t15;
			t49 = t10.*t16;
			t51 = t10.*t20;
			t54 = t11.*t22;
			t50 = -t24;
			t52 = -t32;
			t53 = -t33;
			t55 = -t43;
			t56 = t3.*t39;
			t57 = t5.*t38;
			t58 = -t46;
			t59 = t12.*t36;
			t60 = t15+t45;
			t61 = t16+t44;
			t62 = t28+t31;
			t63 = t29+t30;
			t70 = t42+t43;
			t72 = t24+t54;
			t64 = qdot1.*t60;
			t65 = qdot1.*t61;
			t66 = qdot3.*t60;
			t67 = qdot3.*t61;
			t68 = t12+t58;
			t69 = t20+t55;
			t71 = t23+t52;
			t73 = t25+t53;
			t74 = t34+t50;
			t83 = t39+t57;
			t84 = t40+t56;
			R = reshape([t68,t21,t84,t19,t3,t14,t60,t41,t70],[3,3]);
			if nargout > 1
				dRdq = reshape([t84,0.0,-t12+t46,t14,0.0,-t19,t70,0.0,t83,t7.*t19,t17,t7.*t14,t3.*t5,t22,t2.*t3,t16.*t22,t38,t12.*t22,t83,t18,t69,0.0,0.0,0.0,t68,t21,t84],[3,3,3]);
			end
			if nargout > 2
				t75 = qdot1.*t68;
				t76 = qdot1.*t69;
				t77 = -t64;
				t78 = -t65;
				t79 = qdot3.*t68;
				t80 = qdot3.*t69;
				t81 = -t66;
				t82 = -t67;
				t85 = -t75;
				t86 = -t79;
				t87 = t51+t78+t81;
				t88 = t59+t77+t82;
				t89 = t49+t76+t86;
				t90 = t48+t80+t85;
				Rdot = reshape([t87,t62,t90,t26+t27,t36,t71,-t76+t79+t16.*t36,t74,t88],[3,3]);
			end
			if nargout > 3
				dRdotdq = reshape([t90,0.0,t65+t66+t20.*t36,t71,0.0,-t26-t27,t88,0.0,t89,t8.*t20+t9.*t19+t7.*t26,t73,t8.*t15+t9.*t14+t7.*t52,t5.*t36+qdot1.*t2.*t3,t35,t2.*t36-qdot1.*t3.*t5,t11.*t19+t16.*t35+qdot1.*t12.*t22,t63,t11.*t14+t12.*t35+qdot1.*t6.*t16,t89,t72,t47+t64+t67,0.0,0.0,0.0,t87,t62,t90],[3,3,3]);
			end
			if nargout > 4
				S = reshape([t21,t3,t41,t4,0.0,t7,0.0,1.0,0.0],[3,3]);
			end
			if nargout > 5
				detS = t22;
			end
			if nargout > 6
				dSdq = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t17,t22,t38,0.0,0.0,0.0,0.0,0.0,0.0,t18,0.0,t21,-t7,0.0,t4,0.0,0.0,0.0],[3,3,3]);
			end
			if nargout > 7
				Sdot = reshape([t62,t36,t74,t37,0.0,t9,0.0,0.0,0.0],[3,3]);
			end
			if nargout > 8
				dSdotdq = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t73,t35,t63,0.0,0.0,0.0,0.0,0.0,0.0,t72,0.0,t62,-t9,0.0,t37,0.0,0.0,0.0],[3,3,3]);
			end
		end
		
		%%
		function [R,dRdq,Rdot,dRdotdq,S,detS,dSdq,Sdot,dSdotdq] = ZXZ(in1,in2)
			q1 = in1(1,:);
			q2 = in1(2,:);
			q3 = in1(3,:);
			qdot1 = in2(1,:);
			qdot2 = in2(2,:);
			qdot3 = in2(3,:);
			t2 = cos(q1);
			t3 = cos(q2);
			t4 = cos(q3);
			t5 = sin(q1);
			t6 = sin(q2);
			t7 = sin(q3);
			t8 = qdot2.*t3;
			t9 = qdot3.*t4;
			t10 = qdot2.*t6;
			t11 = qdot3.*t7;
			t12 = t2.*t4;
			t13 = t3.*t4;
			t14 = t2.*t6;
			t15 = t2.*t7;
			t16 = t4.*t5;
			t17 = t3.*t7;
			t18 = t4.*t6;
			t19 = t5.*t6;
			t20 = t5.*t7;
			t21 = t6.*t7;
			t22 = -t6;
			t23 = -t7;
			t24 = t2.*t8;
			t25 = t4.*t8;
			t26 = t3.*t9;
			t27 = qdot1.*t14;
			t28 = t5.*t8;
			t29 = t7.*t8;
			t30 = t4.*t10;
			t31 = t3.*t11;
			t32 = t6.*t9;
			t33 = qdot1.*t19;
			t36 = -t8;
			t37 = -t9;
			t38 = -t10;
			t39 = -t15;
			t40 = -t16;
			t41 = -t20;
			t42 = -t21;
			t43 = t3.*t12;
			t44 = t3.*t15;
			t45 = t5.*t13;
			t46 = t5.*t17;
			t47 = t10.*t12;
			t48 = t10.*t15;
			t49 = t10.*t16;
			t50 = t10.*t20;
			t51 = t8.*t23;
			t54 = t9.*t22;
			t55 = t10.*t23;
			t56 = t11.*t22;
			t52 = -t30;
			t53 = -t31;
			t57 = -t43;
			t58 = t3.*t39;
			t59 = -t45;
			t60 = -t46;
			t61 = t12.*t38;
			t62 = t15+t45;
			t63 = t16+t44;
			t64 = t27+t28;
			t65 = t29+t32;
			t72 = t41+t43;
			t73 = t25+t56;
			t74 = t26+t55;
			t85 = t51+t54;
			t66 = qdot1.*t62;
			t67 = qdot1.*t63;
			t68 = qdot3.*t62;
			t69 = qdot3.*t63;
			t70 = t12+t60;
			t71 = t20+t57;
			t83 = t39+t59;
			R = reshape([t70,t63,t21,t83,t72,t18,t19,-t14,t3],[3,3]);
			if nargout > 1
				t84 = t40+t58;
				dRdq = reshape([t84,t70,0.0,t71,t83,0.0,t14,t19,0.0,t7.*t19,t14.*t23,t17,t6.*t16,t12.*t22,t13,t3.*t5,-t2.*t3,t22,t83,t72,t18,-t12+t46,t84,t42,0.0,0.0,0.0],[3,3,3]);
			end
			if nargout > 2
				t86 = t52+t53;
				t75 = qdot1.*t70;
				t76 = qdot1.*t71;
				t77 = -t66;
				t78 = -t67;
				t79 = qdot3.*t70;
				t80 = qdot3.*t71;
				t81 = -t68;
				t82 = -t69;
				t87 = -t75;
				t88 = -t79;
				t89 = t50+t78+t81;
				t90 = t61+t77+t82;
				t91 = t49+t76+t88;
				Rdot = reshape([t89,t75-t80+t15.*t38,t65,t91,t90,t73,t64,-t24+t33,t38],[3,3]);
			end
			if nargout > 3
				t92 = t48+t80+t87;
				dRdotdq = reshape([t92,t89,0.0,t47+t66+t69,t91,0.0,t24-t33,t64,0.0,t8.*t20+t9.*t19+t7.*t27,t7.*t33+t14.*t37+t15.*t36,t74,t8.*t16-t11.*t19+qdot1.*t6.*t12,t11.*t14+t12.*t36+qdot1.*t6.*t16,t86,t5.*t38+qdot1.*t2.*t3,t2.*t10+qdot1.*t3.*t5,t36,t91,t90,t73,t67+t68+t20.*t38,t92,t85,0.0,0.0,0.0],[3,3,3]);
			end
			if nargout > 4
				S = reshape([t21,t18,t3,t4,t23,0.0,0.0,0.0,1.0],[3,3]);
			end
			if nargout > 5
				detS = t22;
			end
			if nargout > 6
				dSdq = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t17,t13,t22,0.0,0.0,0.0,0.0,0.0,0.0,t18,t42,0.0,t23,-t4,0.0,0.0,0.0,0.0],[3,3,3]);
			end
			if nargout > 7
				Sdot = reshape([t65,t73,t38,-t11,t37,0.0,0.0,0.0,0.0],[3,3]);
			end
			if nargout > 8
				dSdotdq = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t74,t86,t36,0.0,0.0,0.0,0.0,0.0,0.0,t73,t85,0.0,t37,t11,0.0,0.0,0.0,0.0],[3,3,3]);
			end
		end
		
		%%
		function [R,dRdq,Rdot,dRdotdq,S,detS,dSdq,Sdot,dSdotdq] = ZYZ(in1,in2)
			q1 = in1(1,:);
			q2 = in1(2,:);
			q3 = in1(3,:);
			qdot1 = in2(1,:);
			qdot2 = in2(2,:);
			qdot3 = in2(3,:);
			t2 = cos(q1);
			t3 = cos(q2);
			t4 = cos(q3);
			t5 = sin(q1);
			t6 = sin(q2);
			t7 = sin(q3);
			t8 = qdot2.*t3;
			t9 = qdot3.*t4;
			t10 = qdot2.*t6;
			t11 = qdot3.*t7;
			t12 = t2.*t4;
			t13 = t3.*t4;
			t14 = t2.*t6;
			t15 = t2.*t7;
			t16 = t4.*t5;
			t17 = t3.*t7;
			t18 = t4.*t6;
			t19 = t5.*t6;
			t20 = t5.*t7;
			t21 = t6.*t7;
			t22 = -t6;
			t23 = t2.*t8;
			t24 = t4.*t8;
			t25 = t3.*t9;
			t26 = qdot1.*t14;
			t27 = t5.*t8;
			t28 = t7.*t8;
			t29 = t4.*t10;
			t30 = t3.*t11;
			t31 = t6.*t9;
			t32 = qdot1.*t19;
			t33 = t7.*t10;
			t34 = t6.*t11;
			t35 = -t8;
			t36 = -t10;
			t37 = -t11;
			t38 = -t13;
			t39 = -t15;
			t40 = -t16;
			t41 = -t18;
			t42 = -t20;
			t43 = t3.*t12;
			t44 = t3.*t15;
			t45 = t5.*t13;
			t46 = t5.*t17;
			t47 = t10.*t12;
			t48 = t10.*t15;
			t49 = t10.*t16;
			t51 = t10.*t20;
			t54 = t11.*t22;
			t50 = -t24;
			t52 = -t32;
			t53 = -t33;
			t55 = -t43;
			t56 = t3.*t39;
			t57 = t5.*t38;
			t58 = -t46;
			t59 = t12.*t36;
			t60 = t15+t45;
			t61 = t16+t44;
			t62 = t28+t31;
			t63 = t29+t30;
			t70 = t42+t43;
			t72 = t24+t54;
			t64 = qdot1.*t60;
			t65 = qdot1.*t61;
			t66 = qdot3.*t60;
			t67 = qdot3.*t61;
			t68 = t12+t58;
			t69 = t20+t55;
			t71 = t23+t52;
			t73 = t25+t53;
			t74 = t34+t50;
			t83 = t39+t57;
			t84 = t40+t56;
			R = reshape([t70,t60,t41,t84,t68,t21,t14,t19,t3],[3,3]);
			if nargout > 1
				dRdq = reshape([t83,t70,0.0,-t12+t46,t84,0.0,-t19,t14,0.0,t12.*t22,t16.*t22,t38,t7.*t14,t7.*t19,t17,t2.*t3,t3.*t5,t22,t84,t68,t21,t69,t83,t18,0.0,0.0,0.0],[3,3,3]);
			end
			if nargout > 2
				t75 = qdot1.*t68;
				t76 = qdot1.*t69;
				t77 = -t64;
				t78 = -t65;
				t79 = qdot3.*t68;
				t80 = qdot3.*t69;
				t81 = -t66;
				t82 = -t67;
				t85 = -t75;
				t86 = -t79;
				t87 = t51+t78+t81;
				t88 = t59+t77+t82;
				t89 = t49+t76+t86;
				t90 = t48+t80+t85;
				Rdot = reshape([t88,-t76+t79+t16.*t36,t74,t90,t87,t62,t71,t26+t27,t36],[3,3]);
			end
			if nargout > 3
				dRdotdq = reshape([t89,t88,0.0,t65+t66+t20.*t36,t90,0.0,-t26-t27,t71,0.0,t11.*t14+t12.*t35+qdot1.*t6.*t16,t11.*t19+t16.*t35+qdot1.*t12.*t22,t63,t8.*t15+t9.*t14+t7.*t52,t8.*t20+t9.*t19+t7.*t26,t73,t2.*t36-qdot1.*t3.*t5,t5.*t36+qdot1.*t2.*t3,t35,t90,t87,t62,t47+t64+t67,t89,t72,0.0,0.0,0.0],[3,3,3]);
			end
			if nargout > 4
				S = reshape([t41,t21,t3,t7,t4,0.0,0.0,0.0,1.0],[3,3]);
			end
			if nargout > 5
				detS = t22;
			end
			if nargout > 6
				dSdq = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t38,t17,t22,0.0,0.0,0.0,0.0,0.0,0.0,t21,t18,0.0,t4,-t7,0.0,0.0,0.0,0.0],[3,3,3]);
			end
			if nargout > 7
				Sdot = reshape([t74,t62,t36,t9,t37,0.0,0.0,0.0,0.0],[3,3]);
			end
			if nargout > 8
				dSdotdq = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t63,t73,t35,0.0,0.0,0.0,0.0,0.0,0.0,t62,t72,0.0,t37,-t9,0.0,0.0,0.0,0.0],[3,3,3]);
			end
		end
		
		%%
		function [R,dRdq,Rdot,dRdotdq,S,detS,dSdq,Sdot,dSdotdq] = XYZ(in1,in2)
			q1 = in1(1,:);
			q2 = in1(2,:);
			q3 = in1(3,:);
			qdot1 = in2(1,:);
			qdot2 = in2(2,:);
			qdot3 = in2(3,:);
			t2 = cos(q1);
			t3 = cos(q2);
			t4 = cos(q3);
			t5 = sin(q1);
			t6 = sin(q2);
			t7 = sin(q3);
			t8 = qdot2.*t3;
			t9 = qdot3.*t4;
			t10 = qdot2.*t6;
			t11 = qdot3.*t7;
			t12 = t2.*t3;
			t13 = t2.*t4;
			t14 = t3.*t4;
			t15 = t3.*t5;
			t16 = t2.*t7;
			t17 = t4.*t5;
			t18 = t3.*t7;
			t19 = t4.*t6;
			t20 = t5.*t7;
			t21 = t6.*t7;
			t22 = qdot1.*t12;
			t23 = t4.*t8;
			t24 = t3.*t9;
			t25 = qdot1.*t15;
			t26 = t2.*t10;
			t27 = t7.*t8;
			t28 = t4.*t10;
			t29 = t3.*t11;
			t30 = t6.*t9;
			t31 = t5.*t10;
			t32 = t7.*t10;
			t33 = t6.*t11;
			t34 = -t10;
			t35 = -t11;
			t36 = -t14;
			t37 = t6.*t20;
			t38 = -t15;
			t39 = -t18;
			t40 = -t19;
			t41 = -t20;
			t42 = t6.*t13;
			t43 = t6.*t16;
			t44 = t6.*t17;
			t45 = t8.*t13;
			t46 = t8.*t16;
			t47 = t8.*t17;
			t49 = -t22;
			t50 = -t23;
			t51 = -t24;
			t52 = -t28;
			t53 = -t29;
			t54 = -t37;
			t55 = -t42;
			t56 = t8.*t41;
			t57 = t27+t30;
			t58 = t28+t29;
			t59 = t16+t44;
			t60 = t17+t43;
			t70 = t41+t42;
			t61 = qdot1.*t59;
			t62 = qdot1.*t60;
			t63 = qdot3.*t59;
			t64 = qdot3.*t60;
			t65 = t31+t49;
			t66 = t33+t50;
			t67 = t32+t51;
			t68 = t13+t54;
			t69 = t20+t55;
			R = reshape([t14,t59,t69,t39,t68,t60,t6,t38,t12],[3,3]);
			if nargout > 1
				dRdq = reshape([0.0,t70,t59,0.0,-t60,t68,0.0,-t12,t38,t40,t5.*t14,-t4.*t12,t21,t7.*t38,t7.*t12,t3,t5.*t6,-t2.*t6,t39,t68,t60,t36,-t59,t70,0.0,0.0,0.0],[3,3,3]);
			end
			if nargout > 2
				t79 = t52+t53;
				t71 = qdot1.*t68;
				t72 = qdot1.*t69;
				t73 = -t61;
				t74 = -t62;
				t75 = qdot3.*t68;
				t76 = qdot3.*t69;
				t77 = -t63;
				t78 = -t64;
				t80 = -t72;
				t81 = -t76;
				t82 = t45+t73+t78;
				t85 = t56+t74+t77;
				t83 = t46+t71+t81;
				t84 = t47+t75+t80;
				Rdot = reshape([t79,t84,-t45+t61+t64,t67,t85,t83,t8,t65,-t25-t26],[3,3]);
			end
			if nargout > 3
				dRdotdq = reshape([0.0,t82,t84,0.0,-t46-t71+t76,t85,0.0,t25+t26,t65,t66,t4.*t22+t15.*t35+t17.*t34,t10.*t13+t11.*t12+qdot1.*t5.*t14,t57,t10.*t20+t9.*t38+t7.*t49,t9.*t12-t7.*t25+t16.*t34,t34,t5.*t8+qdot1.*t2.*t6,-t2.*t8+qdot1.*t5.*t6,t67,t85,t83,t58,-t47+t72-t75,t82,0.0,0.0,0.0],[3,3,3]);
			end
			if nargout > 4
				S = reshape([t14,t39,t6,t7,t4,0.0,0.0,0.0,1.0],[3,3]);
			end
			if nargout > 5
				detS = t3;
			end
			if nargout > 6
				dSdq = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t40,t21,t3,0.0,0.0,0.0,0.0,0.0,0.0,t39,t36,0.0,t4,-t7,0.0,0.0,0.0,0.0],[3,3,3]);
			end
			if nargout > 7
				Sdot = reshape([t79,t67,t8,t9,t35,0.0,0.0,0.0,0.0],[3,3]);
			end
			if nargout > 8
				dSdotdq = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t66,t57,t34,0.0,0.0,0.0,0.0,0.0,0.0,t67,t58,0.0,t35,-t9,0.0,0.0,0.0,0.0],[3,3,3]);
			end
		end
		
		%%
		function [R,dRdq,Rdot,dRdotdq,S,detS,dSdq,Sdot,dSdotdq] = XZY(in1,in2)
			q1 = in1(1,:);
			q2 = in1(2,:);
			q3 = in1(3,:);
			qdot1 = in2(1,:);
			qdot2 = in2(2,:);
			qdot3 = in2(3,:);
			t2 = cos(q1);
			t3 = cos(q2);
			t4 = cos(q3);
			t5 = sin(q1);
			t6 = sin(q2);
			t7 = sin(q3);
			t8 = qdot2.*t3;
			t9 = qdot3.*t4;
			t10 = qdot2.*t6;
			t11 = qdot3.*t7;
			t12 = t2.*t3;
			t13 = t2.*t4;
			t14 = t3.*t4;
			t15 = t3.*t5;
			t16 = t2.*t7;
			t17 = t4.*t5;
			t18 = t3.*t7;
			t19 = t4.*t6;
			t20 = t5.*t7;
			t21 = t6.*t7;
			t22 = -t3;
			t23 = -t6;
			t24 = -t7;
			t25 = qdot1.*t12;
			t26 = t4.*t8;
			t27 = t3.*t9;
			t28 = qdot1.*t15;
			t29 = t2.*t10;
			t31 = t4.*t10;
			t34 = t5.*t10;
			t35 = t7.*t10;
			t36 = t6.*t11;
			t37 = -t8;
			t38 = -t9;
			t39 = -t13;
			t40 = t6.*t20;
			t41 = -t16;
			t42 = -t17;
			t43 = -t18;
			t44 = -t19;
			t45 = -t21;
			t46 = t6.*t13;
			t47 = t6.*t16;
			t48 = t6.*t17;
			t49 = t8.*t13;
			t50 = t8.*t16;
			t51 = t8.*t17;
			t52 = t8.*t20;
			t54 = t9.*t22;
			t57 = t8.*t24;
			t59 = t11.*t22;
			t60 = t9.*t23;
			t61 = t10.*t24;
			t62 = t20.*t23;
			t63 = t16.*t23;
			t64 = t17.*t23;
			t53 = -t26;
			t55 = -t28;
			t56 = -t29;
			t58 = -t31;
			t65 = t20.*t37;
			t66 = t13+t40;
			t67 = t20+t46;
			t72 = t27+t61;
			t74 = t35+t54;
			t75 = t16+t64;
			t76 = t17+t63;
			t77 = t41+t48;
			t78 = t42+t47;
			R = reshape([t14,t67,t77,t23,t12,t15,t18,t78,t66],[3,3]);
			if nargout > 1
				t86 = t57+t60;
				t88 = t39+t62;
				dRdq = reshape([0.0,t75,t67,0.0,-t15,t12,0.0,t88,t78,t44,t4.*t12,t5.*t14,t22,t2.*t23,t5.*t23,t45,t7.*t12,t7.*t15,t43,t76,t88,0.0,0.0,0.0,t14,t67,t77],[3,3,3]);
			end
			if nargout > 2
				t68 = qdot1.*t66;
				t69 = qdot1.*t67;
				t70 = qdot3.*t66;
				t71 = qdot3.*t67;
				t73 = t36+t53;
				t79 = qdot1.*t75;
				t80 = qdot1.*t76;
				t82 = qdot3.*t75;
				t83 = qdot3.*t76;
				t85 = t55+t56;
				t87 = t58+t59;
				t81 = -t68;
				t84 = -t70;
				t91 = t49+t79+t83;
				t92 = t65+t80+t82;
				t89 = t51+t69+t84;
				t90 = t50+t71+t81;
				Rdot = reshape([t87,t91,t89,t37,t85,t25-t34,t72,t90,t52-t80-t82],[3,3]);
			end
			if nargout > 3
				dRdotdq = reshape([0.0,-t69+t70+t17.*t37,t91,0.0,-t25+t34,t85,0.0,t92,t90,t73,-t11.*t12+t10.*t39-qdot1.*t5.*t14,-t11.*t15+t4.*t25+t10.*t42,t10,t2.*t37+qdot1.*t5.*t6,t5.*t37+qdot1.*t2.*t23,t86,t9.*t12+t10.*t41+t24.*t28,t9.*t15-t10.*t20+t7.*t25,t74,t68-t71+t16.*t37,t92,0.0,0.0,0.0,t87,t91,t89],[3,3,3]);
			end
			if nargout > 4
				S = reshape([t14,t23,t18,t24,0.0,t4,0.0,1.0,0.0],[3,3]);
			end
			if nargout > 5
				detS = t22;
			end
			if nargout > 6
				dSdq = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t44,t22,t45,0.0,0.0,0.0,0.0,0.0,0.0,t43,0.0,t14,-t4,0.0,t24,0.0,0.0,0.0],[3,3,3]);
			end
			if nargout > 7
				Sdot = reshape([t87,t37,t72,t38,0.0,-t11,0.0,0.0,0.0],[3,3]);
			end
			if nargout > 8
				dSdotdq = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t73,t10,t86,0.0,0.0,0.0,0.0,0.0,0.0,t74,0.0,t87,t11,0.0,t38,0.0,0.0,0.0],[3,3,3]);
			end
		end
		
		%%
		function [R,dRdq,Rdot,dRdotdq,S,detS,dSdq,Sdot,dSdotdq] = YZX(in1,in2)
			q1 = in1(1,:);
			q2 = in1(2,:);
			q3 = in1(3,:);
			qdot1 = in2(1,:);
			qdot2 = in2(2,:);
			qdot3 = in2(3,:);
			t2 = cos(q1);
			t3 = cos(q2);
			t4 = cos(q3);
			t5 = sin(q1);
			t6 = sin(q2);
			t7 = sin(q3);
			t8 = qdot2.*t3;
			t9 = qdot3.*t4;
			t10 = qdot2.*t6;
			t11 = qdot3.*t7;
			t12 = t2.*t3;
			t13 = t2.*t4;
			t14 = t3.*t4;
			t15 = t3.*t5;
			t16 = t2.*t7;
			t17 = t4.*t5;
			t18 = t3.*t7;
			t19 = t4.*t6;
			t20 = t5.*t7;
			t21 = t6.*t7;
			t22 = qdot1.*t12;
			t23 = t4.*t8;
			t24 = t3.*t9;
			t25 = qdot1.*t15;
			t26 = t2.*t10;
			t27 = t7.*t8;
			t28 = t4.*t10;
			t29 = t3.*t11;
			t30 = t6.*t9;
			t31 = t5.*t10;
			t32 = t7.*t10;
			t33 = t6.*t11;
			t34 = -t10;
			t35 = -t11;
			t36 = -t14;
			t37 = t6.*t20;
			t38 = -t15;
			t39 = -t18;
			t40 = -t19;
			t41 = -t20;
			t42 = t6.*t13;
			t43 = t6.*t16;
			t44 = t6.*t17;
			t45 = t8.*t13;
			t46 = t8.*t16;
			t47 = t8.*t17;
			t49 = -t22;
			t50 = -t23;
			t51 = -t24;
			t52 = -t28;
			t53 = -t29;
			t54 = -t37;
			t55 = -t42;
			t56 = t8.*t41;
			t57 = t27+t30;
			t58 = t28+t29;
			t59 = t16+t44;
			t60 = t17+t43;
			t70 = t41+t42;
			t61 = qdot1.*t59;
			t62 = qdot1.*t60;
			t63 = qdot3.*t59;
			t64 = qdot3.*t60;
			t65 = t31+t49;
			t66 = t33+t50;
			t67 = t32+t51;
			t68 = t13+t54;
			t69 = t20+t55;
			R = reshape([t12,t6,t38,t69,t14,t59,t60,t39,t68],[3,3]);
			if nargout > 1
				dRdq = reshape([t38,0.0,-t12,t59,0.0,t70,t68,0.0,-t60,-t2.*t6,t3,t5.*t6,-t4.*t12,t40,t5.*t14,t7.*t12,t21,t7.*t38,0.0,0.0,0.0,t60,t39,t68,t70,t36,-t59],[3,3,3]);
			end
			if nargout > 2
				t79 = t52+t53;
				t71 = qdot1.*t68;
				t72 = qdot1.*t69;
				t73 = -t61;
				t74 = -t62;
				t75 = qdot3.*t68;
				t76 = qdot3.*t69;
				t77 = -t63;
				t78 = -t64;
				t80 = -t72;
				t81 = -t76;
				t82 = t45+t73+t78;
				t85 = t56+t74+t77;
				t83 = t46+t71+t81;
				t84 = t47+t75+t80;
				Rdot = reshape([-t25-t26,t8,t65,-t45+t61+t64,t79,t84,t83,t67,t85],[3,3]);
			end
			if nargout > 3
				dRdotdq = reshape([t65,0.0,t25+t26,t84,0.0,t82,t85,0.0,-t46-t71+t76,-t2.*t8+qdot1.*t5.*t6,t34,t5.*t8+qdot1.*t2.*t6,t10.*t13+t11.*t12+qdot1.*t5.*t14,t66,t4.*t22+t15.*t35+t17.*t34,t9.*t12-t7.*t25+t16.*t34,t57,t10.*t20+t9.*t38+t7.*t49,0.0,0.0,0.0,t83,t67,t85,t82,t58,-t47+t72-t75],[3,3,3]);
			end
			if nargout > 4
				S = reshape([t6,t14,t39,0.0,t7,t4,1.0,0.0,0.0],[3,3]);
			end
			if nargout > 5
				detS = t3;
			end
			if nargout > 6
				dSdq = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t3,t40,t21,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t39,t36,0.0,t4,-t7,0.0,0.0,0.0],[3,3,3]);
			end
			if nargout > 7
				Sdot = reshape([t8,t79,t67,0.0,t9,t35,0.0,0.0,0.0],[3,3]);
			end
			if nargout > 8
				dSdotdq = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t34,t66,t57,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t67,t58,0.0,t35,-t9,0.0,0.0,0.0],[3,3,3]);
			end
		end
		
		%%
		function [R,dRdq,Rdot,dRdotdq,S,detS,dSdq,Sdot,dSdotdq] = YXZ(in1,in2)
			q1 = in1(1,:);
			q2 = in1(2,:);
			q3 = in1(3,:);
			qdot1 = in2(1,:);
			qdot2 = in2(2,:);
			qdot3 = in2(3,:);
			t2 = cos(q1);
			t3 = cos(q2);
			t4 = cos(q3);
			t5 = sin(q1);
			t6 = sin(q2);
			t7 = sin(q3);
			t8 = qdot2.*t3;
			t9 = qdot3.*t4;
			t10 = qdot2.*t6;
			t11 = qdot3.*t7;
			t12 = t2.*t3;
			t13 = t2.*t4;
			t14 = t3.*t4;
			t15 = t3.*t5;
			t16 = t2.*t7;
			t17 = t4.*t5;
			t18 = t3.*t7;
			t19 = t4.*t6;
			t20 = t5.*t7;
			t21 = t6.*t7;
			t22 = -t3;
			t23 = -t6;
			t24 = -t7;
			t25 = qdot1.*t12;
			t26 = t4.*t8;
			t27 = t3.*t9;
			t28 = qdot1.*t15;
			t29 = t2.*t10;
			t31 = t4.*t10;
			t34 = t5.*t10;
			t35 = t7.*t10;
			t36 = t6.*t11;
			t37 = -t8;
			t38 = -t9;
			t39 = -t13;
			t40 = t6.*t20;
			t41 = -t16;
			t42 = -t17;
			t43 = -t18;
			t44 = -t19;
			t45 = -t21;
			t46 = t6.*t13;
			t47 = t6.*t16;
			t48 = t6.*t17;
			t49 = t8.*t13;
			t50 = t8.*t16;
			t51 = t8.*t17;
			t52 = t8.*t20;
			t54 = t9.*t22;
			t57 = t8.*t24;
			t59 = t11.*t22;
			t60 = t9.*t23;
			t61 = t10.*t24;
			t62 = t20.*t23;
			t63 = t16.*t23;
			t64 = t17.*t23;
			t53 = -t26;
			t55 = -t28;
			t56 = -t29;
			t58 = -t31;
			t65 = t20.*t37;
			t66 = t13+t40;
			t67 = t20+t46;
			t72 = t27+t61;
			t74 = t35+t54;
			t75 = t16+t64;
			t76 = t17+t63;
			t77 = t41+t48;
			t78 = t42+t47;
			R = reshape([t66,t18,t78,t77,t14,t67,t15,t23,t12],[3,3]);
			if nargout > 1
				t86 = t57+t60;
				t88 = t39+t62;
				dRdq = reshape([t78,0.0,t88,t67,0.0,t75,t12,0.0,-t15,t7.*t15,t45,t7.*t12,t5.*t14,t44,t4.*t12,t5.*t23,t22,t2.*t23,t77,t14,t67,t88,t43,t76,0.0,0.0,0.0],[3,3,3]);
			end
			if nargout > 2
				t68 = qdot1.*t66;
				t69 = qdot1.*t67;
				t70 = qdot3.*t66;
				t71 = qdot3.*t67;
				t73 = t36+t53;
				t79 = qdot1.*t75;
				t80 = qdot1.*t76;
				t82 = qdot3.*t75;
				t83 = qdot3.*t76;
				t85 = t55+t56;
				t87 = t58+t59;
				t81 = -t68;
				t84 = -t70;
				t91 = t49+t79+t83;
				t92 = t65+t80+t82;
				t89 = t51+t69+t84;
				t90 = t50+t71+t81;
				Rdot = reshape([t52-t80-t82,t72,t90,t89,t87,t91,t25-t34,t37,t85],[3,3]);
			end
			if nargout > 3
				dRdotdq = reshape([t90,0.0,t92,t91,0.0,-t69+t70+t17.*t37,t85,0.0,-t25+t34,t9.*t15-t10.*t20+t7.*t25,t86,t9.*t12+t10.*t41+t24.*t28,-t11.*t15+t4.*t25+t10.*t42,t73,-t11.*t12+t10.*t39-qdot1.*t5.*t14,t5.*t37+qdot1.*t2.*t23,t10,t2.*t37+qdot1.*t5.*t6,t89,t87,t91,t92,t74,t68-t71+t16.*t37,0.0,0.0,0.0],[3,3,3]);
			end
			if nargout > 4
				S = reshape([t18,t14,t23,t4,t24,0.0,0.0,0.0,1.0],[3,3]);
			end
			if nargout > 5
				detS = t22;
			end
			if nargout > 6
				dSdq = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t45,t44,t22,0.0,0.0,0.0,0.0,0.0,0.0,t14,t43,0.0,t24,-t4,0.0,0.0,0.0,0.0],[3,3,3]);
			end
			if nargout > 7
				Sdot = reshape([t72,t87,t37,-t11,t38,0.0,0.0,0.0,0.0],[3,3]);
			end
			if nargout > 8
				dSdotdq = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t86,t73,t10,0.0,0.0,0.0,0.0,0.0,0.0,t87,t74,0.0,t38,t11,0.0,0.0,0.0,0.0],[3,3,3]);
			end
		end
		
		%%
		function [R,dRdq,Rdot,dRdotdq,S,detS,dSdq,Sdot,dSdotdq] = ZXY(in1,in2)
			q1 = in1(1,:);
			q2 = in1(2,:);
			q3 = in1(3,:);
			qdot1 = in2(1,:);
			qdot2 = in2(2,:);
			qdot3 = in2(3,:);
			t2 = cos(q1);
			t3 = cos(q2);
			t4 = cos(q3);
			t5 = sin(q1);
			t6 = sin(q2);
			t7 = sin(q3);
			t8 = qdot2.*t3;
			t9 = qdot3.*t4;
			t10 = qdot2.*t6;
			t11 = qdot3.*t7;
			t12 = t2.*t3;
			t13 = t2.*t4;
			t14 = t3.*t4;
			t15 = t3.*t5;
			t16 = t2.*t7;
			t17 = t4.*t5;
			t18 = t3.*t7;
			t19 = t4.*t6;
			t20 = t5.*t7;
			t21 = t6.*t7;
			t22 = qdot1.*t12;
			t23 = t4.*t8;
			t24 = t3.*t9;
			t25 = qdot1.*t15;
			t26 = t2.*t10;
			t27 = t7.*t8;
			t28 = t4.*t10;
			t29 = t3.*t11;
			t30 = t6.*t9;
			t31 = t5.*t10;
			t32 = t7.*t10;
			t33 = t6.*t11;
			t34 = -t10;
			t35 = -t11;
			t36 = -t14;
			t37 = t6.*t20;
			t38 = -t15;
			t39 = -t18;
			t40 = -t19;
			t41 = -t20;
			t42 = t6.*t13;
			t43 = t6.*t16;
			t44 = t6.*t17;
			t45 = t8.*t13;
			t46 = t8.*t16;
			t47 = t8.*t17;
			t49 = -t22;
			t50 = -t23;
			t51 = -t24;
			t52 = -t28;
			t53 = -t29;
			t54 = -t37;
			t55 = -t42;
			t56 = t8.*t41;
			t57 = t27+t30;
			t58 = t28+t29;
			t59 = t16+t44;
			t60 = t17+t43;
			t70 = t41+t42;
			t61 = qdot1.*t59;
			t62 = qdot1.*t60;
			t63 = qdot3.*t59;
			t64 = qdot3.*t60;
			t65 = t31+t49;
			t66 = t33+t50;
			t67 = t32+t51;
			t68 = t13+t54;
			t69 = t20+t55;
			R = reshape([t68,t60,t39,t38,t12,t6,t59,t69,t14],[3,3]);
			if nargout > 1
				dRdq = reshape([-t60,t68,0.0,-t12,t38,0.0,t70,t59,0.0,t7.*t38,t7.*t12,t21,t5.*t6,-t2.*t6,t3,t5.*t14,-t4.*t12,t40,-t59,t70,t36,0.0,0.0,0.0,t68,t60,t39],[3,3,3]);
			end
			if nargout > 2
				t79 = t52+t53;
				t71 = qdot1.*t68;
				t72 = qdot1.*t69;
				t73 = -t61;
				t74 = -t62;
				t75 = qdot3.*t68;
				t76 = qdot3.*t69;
				t77 = -t63;
				t78 = -t64;
				t80 = -t72;
				t81 = -t76;
				t82 = t45+t73+t78;
				t85 = t56+t74+t77;
				t83 = t46+t71+t81;
				t84 = t47+t75+t80;
				Rdot = reshape([t85,t83,t67,t65,-t25-t26,t8,t84,-t45+t61+t64,t79],[3,3]);
			end
			if nargout > 3
				dRdotdq = reshape([-t46-t71+t76,t85,0.0,t25+t26,t65,0.0,t82,t84,0.0,t10.*t20+t9.*t38+t7.*t49,t9.*t12-t7.*t25+t16.*t34,t57,t5.*t8+qdot1.*t2.*t6,-t2.*t8+qdot1.*t5.*t6,t34,t4.*t22+t15.*t35+t17.*t34,t10.*t13+t11.*t12+qdot1.*t5.*t14,t66,-t47+t72-t75,t82,t58,0.0,0.0,0.0,t85,t83,t67],[3,3,3]);
			end
			if nargout > 4
				S = reshape([t39,t6,t14,t4,0.0,t7,0.0,1.0,0.0],[3,3]);
			end
			if nargout > 5
				detS = t3;
			end
			if nargout > 6
				dSdq = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t21,t3,t40,0.0,0.0,0.0,0.0,0.0,0.0,t36,0.0,t39,-t7,0.0,t4,0.0,0.0,0.0],[3,3,3]);
			end
			if nargout > 7
				Sdot = reshape([t67,t8,t79,t35,0.0,t9,0.0,0.0,0.0],[3,3]);
			end
			if nargout > 8
				dSdotdq = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t57,t34,t66,0.0,0.0,0.0,0.0,0.0,0.0,t58,0.0,t67,-t9,0.0,t35,0.0,0.0,0.0],[3,3,3]);
			end
		end
		
		%%
		function [R,dRdq,Rdot,dRdotdq,S,detS,dSdq,Sdot,dSdotdq] = ZYX(in1,in2)
			q1 = in1(1,:);
			q2 = in1(2,:);
			q3 = in1(3,:);
			qdot1 = in2(1,:);
			qdot2 = in2(2,:);
			qdot3 = in2(3,:);
			t2 = cos(q1);
			t3 = cos(q2);
			t4 = cos(q3);
			t5 = sin(q1);
			t6 = sin(q2);
			t7 = sin(q3);
			t8 = qdot2.*t3;
			t9 = qdot3.*t4;
			t10 = qdot2.*t6;
			t11 = qdot3.*t7;
			t12 = t2.*t3;
			t13 = t2.*t4;
			t14 = t3.*t4;
			t15 = t3.*t5;
			t16 = t2.*t7;
			t17 = t4.*t5;
			t18 = t3.*t7;
			t19 = t4.*t6;
			t20 = t5.*t7;
			t21 = t6.*t7;
			t22 = -t3;
			t23 = -t6;
			t24 = -t7;
			t25 = qdot1.*t12;
			t26 = t4.*t8;
			t27 = t3.*t9;
			t28 = qdot1.*t15;
			t29 = t2.*t10;
			t31 = t4.*t10;
			t34 = t5.*t10;
			t35 = t7.*t10;
			t36 = t6.*t11;
			t37 = -t8;
			t38 = -t9;
			t39 = -t13;
			t40 = t6.*t20;
			t41 = -t16;
			t42 = -t17;
			t43 = -t18;
			t44 = -t19;
			t45 = -t21;
			t46 = t6.*t13;
			t47 = t6.*t16;
			t48 = t6.*t17;
			t49 = t8.*t13;
			t50 = t8.*t16;
			t51 = t8.*t17;
			t52 = t8.*t20;
			t54 = t9.*t22;
			t57 = t8.*t24;
			t59 = t11.*t22;
			t60 = t9.*t23;
			t61 = t10.*t24;
			t62 = t20.*t23;
			t63 = t16.*t23;
			t64 = t17.*t23;
			t53 = -t26;
			t55 = -t28;
			t56 = -t29;
			t58 = -t31;
			t65 = t20.*t37;
			t66 = t13+t40;
			t67 = t20+t46;
			t72 = t27+t61;
			t74 = t35+t54;
			t75 = t16+t64;
			t76 = t17+t63;
			t77 = t41+t48;
			t78 = t42+t47;
			R = reshape([t12,t15,t23,t78,t66,t18,t67,t77,t14],[3,3]);
			if nargout > 1
				t86 = t57+t60;
				t88 = t39+t62;
				dRdq = reshape([-t15,t12,0.0,t88,t78,0.0,t75,t67,0.0,t2.*t23,t5.*t23,t22,t7.*t12,t7.*t15,t45,t4.*t12,t5.*t14,t44,0.0,0.0,0.0,t67,t77,t14,t76,t88,t43],[3,3,3]);
			end
			if nargout > 2
				t68 = qdot1.*t66;
				t69 = qdot1.*t67;
				t70 = qdot3.*t66;
				t71 = qdot3.*t67;
				t73 = t36+t53;
				t79 = qdot1.*t75;
				t80 = qdot1.*t76;
				t82 = qdot3.*t75;
				t83 = qdot3.*t76;
				t85 = t55+t56;
				t87 = t58+t59;
				t81 = -t68;
				t84 = -t70;
				t91 = t49+t79+t83;
				t92 = t65+t80+t82;
				t89 = t51+t69+t84;
				t90 = t50+t71+t81;
				Rdot = reshape([t85,t25-t34,t37,t90,t52-t80-t82,t72,t91,t89,t87],[3,3]);
			end
			if nargout > 3
				dRdotdq = reshape([-t25+t34,t85,0.0,t92,t90,0.0,-t69+t70+t17.*t37,t91,0.0,t2.*t37+qdot1.*t5.*t6,t5.*t37+qdot1.*t2.*t23,t10,t9.*t12+t10.*t41+t24.*t28,t9.*t15-t10.*t20+t7.*t25,t86,-t11.*t12+t10.*t39-qdot1.*t5.*t14,-t11.*t15+t4.*t25+t10.*t42,t73,0.0,0.0,0.0,t91,t89,t87,t68-t71+t16.*t37,t92,t74],[3,3,3]);
			end
			if nargout > 4
				S = reshape([t23,t18,t14,0.0,t4,t24,1.0,0.0,0.0],[3,3]);
			end
			if nargout > 5
				detS = t22;
			end
			if nargout > 6
				dSdq = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t22,t45,t44,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t14,t43,0.0,t24,-t4,0.0,0.0,0.0],[3,3,3]);
			end
			if nargout > 7
				Sdot = reshape([t37,t72,t87,0.0,-t11,t38,0.0,0.0,0.0],[3,3]);
			end
			if nargout > 8
				dSdotdq = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t10,t86,t73,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t87,t74,0.0,t38,t11,0.0,0.0,0.0],[3,3,3]);
			end
		end
		
		%%
		function q = XYXinv(R)
			r11 = R(1,1);
			if -1.0 < r11 && r11 < 1.0
				q(1,1) = atan2(R(2,1),-R(3,1));
				q(2,1) = acos(r11);
				q(3,1) = atan2(R(1,2),R(1,3));
			else
				% Gimbal lock
				q = nan(3,1);
			end
		end
		
		%%
		function q = XZXinv(R)
			r11 = R(1,1);
			if -1.0 < r11 && r11 < 1.0
				q(1,1) = atan2(R(3,1),R(2,1));
				q(2,1) = acos(r11);
				q(3,1) = atan2(R(1,3),-R(1,2));
			else
				% Gimbal lock
				q = nan(3,1);
			end
		end
		
		%%
		function q = YXYinv(R)
			r22 = R(2,2);
			if -1.0 < r22 && r22 < 1.0
				q(1,1) = atan2(R(1,2),R(3,2));
				q(2,1) = acos(r22);
				q(3,1) = atan2(R(2,1),-R(2,3));
			else
				% Gimbal lock
				q = nan(3,1);
			end
		end
		
		%%
		function q = YZYinv(R)
			r22 = R(2,2);
			if -1.0 < r22 && r22 < 1.0
				q(1,1) = atan2(R(3,2),-R(1,2));
				q(2,1) = acos(r22);
				q(3,1) = atan2(R(2,3),R(2,1));
			else
				% Gimbal lock
				q = nan(3,1);
			end
		end
		
		%%
		function q = ZXZinv(R)
			r33 = R(3,3);
			if -1.0 < r33 && r33 < 1.0
				q(1,1) = atan2(R(1,3),-R(2,3));
				q(2,1) = acos(r33);
				q(3,1) = atan2(R(3,1),R(3,2));
			else
				% Gimbal lock
				q = nan(3,1);
			end
		end
		
		%%
		function q = ZYZinv(R)
			r33 = R(3,3);
			if -1.0 < r33 && r33 < 1.0
				q(1,1) = atan2(R(2,3),R(1,3));
				q(2,1) = acos(r33);
				q(3,1) = atan2(R(3,2),-R(3,1));
			else
				% Gimbal lock
				q = nan(3,1);
			end
		end
		
		%%
		function q = XYZinv(R)
			r13 = R(1,3);
			if -1.0 < r13 && r13 < 1.0
				q(1,1) = atan2(-R(2,3),R(3,3));
				q(2,1) = asin(r13);
				q(3,1) = atan2(-R(1,2),R(1,1));
			else
				% Gimbal lock
				q = nan(3,1);
			end
		end
		
		%%
		function q = XZYinv(R)
			r12 = R(1,2);
			if -1.0 < r12 && r12 < 1.0
				q(1,1) = atan2(R(3,2),R(2,2));
				q(2,1) = asin(-r12);
				q(3,1) = atan2(R(1,3),R(1,1));
			else
				% Gimbal lock
				q = nan(3,1);
			end
		end
		
		%%
		function q = YXZinv(R)
			r23 = R(2,3);
			if -1.0 < r23 && r23 < 1.0
				q(1,1) = atan2(R(1,3),R(3,3));
				q(2,1) = asin(-r23);
				q(3,1) = atan2(R(2,1),R(2,2));
			else
				% Gimbal lock
				q = nan(3,1);
			end
		end
		
		%%
		function q = YZXinv(R)
			r21 = R(2,1);
			if -1.0 < r21 && r21 < 1.0
				q(1,1) = atan2(-R(3,1),R(1,1));
				q(2,1) = asin(r21);
				q(3,1) = atan2(-R(2,3),R(2,2));
			else
				% Gimbal lock
				q = nan(3,1);
			end
		end
		
		%%
		function q = ZXYinv(R)
			r32 = R(3,2);
			if -1.0 < r32 && r32 < 1.0
				q(1,1) = atan2(-R(1,2),R(2,2));
				q(2,1) = asin(r32);
				q(3,1) = atan2(-R(3,1),R(3,3));
			else
				% Gimbal lock
				q = nan(3,1);
			end
		end
		
		%%
		function q = ZYXinv(R)
			r31 = R(3,1);
			if -1.0 < r31 && r31 < 1.0
				q(1,1) = atan2(R(2,1),R(1,1));
				q(2,1) = asin(-r31);
				q(3,1) = atan2(R(3,2),R(3,3));
			else
				% Gimbal lock
				q = nan(3,1);
			end
		end
	end
end
