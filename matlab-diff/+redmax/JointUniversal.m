classdef JointUniversal < redmax.Joint
	%JointUniversal A universal joint allowing rotation about X and Y
	
	%%
	properties
	end
	
	%%
	methods
		%%
		function this = JointUniversal(parent,body)
			this = this@redmax.Joint(parent,body,2);
		end
	end
	
	%%
	methods (Access = protected)
		%%
		function update_(this,deriv)
			q = this.q;
			qdot = this.qdot;
			
			[R,dRdq,Rdot,dRdotdq,S,dSdq,Sdot,dSdotdq] = redmax.JointUniversal.XY(q,qdot);
			
			this.Q(1:3,1:3) = R;
			this.A = se3.Ad(this.Q);
			
			this.Adot(1:3,1:3) = Rdot;
			this.Adot(4:6,4:6) = Rdot;
			this.S(1:3,1:2) = S;
			this.Sdot(1:3,1:2) = Sdot;
			
			if deriv
				for i = 1 : 2
					this.dAdq(1:3,1:3,i) = dRdq(:,:,i);
					this.dAdq(4:6,4:6,i) = dRdq(:,:,i);
					this.dAdotdq(1:3,1:3,i) = dRdotdq(:,:,i);
					this.dAdotdq(4:6,4:6,i) = dRdotdq(:,:,i);
					this.dSdq(1:3,1:2,i) = dSdq(:,:,i);
					this.dSdotdq(1:3,1:2,i) = dSdotdq(:,:,i);
				end
			end
		end
		
		%%
		function draw_(this)
			s = 2*min(this.body.getAxisSize());
			p = this.E_wj(1:3,4);
			x = this.E_wj(1:3,1);
			y = this.E_wj(1:3,2);
			px0 = p - s*x;
			px1 = p + s*x;
			py0 = p - s*y;
			py1 = p + s*y;
			ps = [px0 px1];
			plot3(ps(1,:),ps(2,:),ps(3,:),'r','LineWidth',2);
			ps = [py0 py1];
			plot3(ps(1,:),ps(2,:),ps(3,:),'g','LineWidth',2);
		end
	end
	
	methods (Static)
		%%
		function codegen()
			clear;
			q = sym('q',[2,1],'real');
			qdot = sym('qdot',[2,1],'real');
			s = sin(q);
			c = cos(q);
			X1 = [1 0 0; 0 c(1) -s(1); 0 s(1) c(1)];
			Y2 = [c(2) 0 s(2); 0 1 0; -s(2) 0 c(2)];			
			R = X1*Y2;
			
			% dRdq
			dRdq = 0*sym('dRdq',[3,3,2],'real');
			for i = 1 : 2
				dRdq(:,:,i) = diff(R,q(i));
			end
			% Rdot
			Rdot = 0*sym('Rdot',[3,3],'real');
			for i = 1 : 2
				Rdot = Rdot + dRdq(:,:,i)*qdot(i);
			end
			Rdot = simplify(Rdot);
			% dRdotdq
			dRdotdq = 0*sym('dRdotdq',[3,3,2],'real');
			for i = 1 : 2
				dRdotdq(:,:,i) = diff(Rdot,q(i));
			end
			% S
			S = 0*sym('S',[3,2],'real');
			for i = 1 : 2
				tmp = simplify(R'*diff(R,q(i)));
				S(1:3,i) = [tmp(3,2); tmp(1,3); tmp(2,1)];
			end
			% dSdq
			dSdq = 0*sym('dSdq',[3,2,2],'real');
			for i = 1 : 2
				dSdq(:,:,i) = diff(S,q(i));
			end
			% Sdot
			Sdot = 0*sym('Sdot',[3,2],'real');
			for i = 1 : 2
				Sdot = Sdot + dSdq(:,:,i)*qdot(i);
			end
			% dSdotdq
			dSdotdq = 0*sym('dSdotdq',[3,2,2],'real');
			for i = 1 : 2
				dSdotdq(:,:,i) = diff(Sdot,q(i));
			end
			% Output
			matlabFunction(R,dRdq,Rdot,dRdotdq,...
				S,dSdq,Sdot,dSdotdq,...
				'File','XY.m','Optimize',true,'Vars',{q,qdot},...
				'Outputs',{'R','dRdq','Rdot','dRdotdq',...
				'S','dSdq','Sdot','dSdotdq'});
		end
	end
	
	methods (Static, Access = private)
		%%
		function [R,dRdq,Rdot,dRdotdq,S,dSdq,Sdot,dSdotdq] = XY(in1,in2)
			q1 = in1(1,:);
			q2 = in1(2,:);
			qdot1 = in2(1,:);
			qdot2 = in2(2,:);
			t2 = cos(q1);
			t3 = cos(q2);
			t4 = sin(q1);
			t5 = sin(q2);
			t6 = qdot1.*t2;
			t7 = qdot2.*t3;
			t8 = qdot1.*t4;
			t9 = qdot2.*t5;
			t10 = t2.*t3;
			t11 = t2.*t5;
			t12 = t3.*t4;
			t13 = t4.*t5;
			t14 = -t5;
			t15 = t3.*t6;
			t16 = t2.*t7;
			t17 = t5.*t6;
			t18 = t3.*t8;
			t19 = t2.*t9;
			t20 = t4.*t7;
			t21 = t5.*t8;
			t22 = t4.*t9;
			t23 = -t7;
			t24 = -t8;
			t25 = -t9;
			t26 = -t10;
			t27 = -t11;
			t28 = -t12;
			R = reshape([t3,t13,t27,0.0,t2,t4,t5,t28,t10],[3,3]);
			if nargout > 1
				dRdq = reshape([0.0,t11,t13,0.0,-t4,t2,0.0,t26,t28,t14,t12,t26,0.0,0.0,0.0,t3,t13,t27],[3,3,2]);
			end
			if nargout > 2
				t29 = -t15;
				t30 = -t16;
				t31 = t17+t20;
				t32 = t18+t19;
				t33 = t22+t29;
				t34 = t21+t30;
				Rdot = reshape([t25,t31,t34,0.0,t24,t6,t7,t33,-t32],[3,3]);
			end
			if nargout > 3
				dRdotdq = reshape([0.0,t16+t8.*t14,t31,0.0,-t6,t24,0.0,t32,t33,t23,t15-t22,t32,0.0,0.0,0.0,t25,t31,t34],[3,3,2]);
			end
			if nargout > 4
				S = reshape([t3,0.0,t5,0.0,1.0,0.0],[3,2]);
			end
			if nargout > 5
				dSdq = reshape([0.0,0.0,0.0,0.0,0.0,0.0,t14,0.0,t3,0.0,0.0,0.0],[3,2,2]);
			end
			if nargout > 6
				Sdot = reshape([t25,0.0,t7,0.0,0.0,0.0],[3,2]);
			end
			if nargout > 7
				dSdotdq = reshape([0.0,0.0,0.0,0.0,0.0,0.0,t23,0.0,t25,0.0,0.0,0.0],[3,2,2]);
			end
		end
	end
end
