classdef ConstraintAttachSpring < redmax.Constraint
	% ConstraintAttachSpring Attaches the two ends of a string to rigid
	% bodies
	
	%%
	properties
		spring
	end
	
	%%
	methods
		%%
		function this = ConstraintAttachSpring(spring)
			this = this@redmax.Constraint(6,0,0,0);
			this.name = [spring.name,'-A'];
			this.spring = spring;
		end
	end
	
	%%
	methods (Access = protected)
		%%
		function [Gm,Gmdot,gm,gmdot,gmddot] = computeJacEqM_(this,Gm,Gmdot,gm,gmdot,gmddot)
			% g = E*r - x
			% G = [R*Gamma, -I]
			% Gdot = [R*[w]*Gamma, 0]
			I = eye(3);
			rows0 = this.idxEM(1:3);
			rows1 = this.idxEM(4:6);
			cols0S = this.spring.nodes(1).idxM;
			cols1S = this.spring.nodes(end).idxM;
			body0 = this.spring.body0;
			body1 = this.spring.body1;
			if isempty(body0)
				E0 = eye(4);
				cols0B = [];
			else
				E0 = body0.E_wi;
				cols0B = body0.idxM;
			end
			if isempty(body1)
				E1 = eye(4);
				cols1B = [];
			else
				E1 = body1.E_wi;
				cols1B = body1.idxM;
			end
			R0 = E0(1:3,1:3);
			R1 = E1(1:3,1:3);
			G0 = se3.Gamma(this.spring.r0);
			G1 = se3.Gamma(this.spring.r1);
			if ~isempty(body0)
				W0 = se3.brac(body0.phi(1:3));
				Gm(rows0,cols0B) = R0*G0;
				Gmdot(rows0,cols0B) = R0*W0*G0;
			end
			if ~isempty(body1)
				W1 = se3.brac(body1.phi(1:3));
				Gm(rows1,cols1B) = R1*G1;
				Gmdot(rows1,cols1B) = R1*W1*G1;
			end
			Gm(rows0,cols0S) = -I;
			Gm(rows1,cols1S) = -I;
			gm0 = E0*[this.spring.r0;1] - [this.spring.nodes(  1).x;1];
			gm1 = E1*[this.spring.r1;1] - [this.spring.nodes(end).x;1];
			gm(rows0) = gm0(1:3,1);
			gm(rows1) = gm1(1:3,1);
		end
	end
end
