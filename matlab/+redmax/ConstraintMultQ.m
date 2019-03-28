classdef ConstraintMultQ < redmax.Constraint
	% ConstraintMultQ Constrains qs to be a multiple of each other
	
	%%
	properties
		jointA % One joint
		jointB % Another joint
		factor % Multiplication factor
	end
	
	%%
	methods
		%%
		function this = ConstraintMultQ(jointA,jointB)
			this = this@redmax.Constraint(0,1,0,0);
			this.name = [jointA.name,'-Q-',jointB.name];
			this.jointA = jointA;
			this.jointB = jointB;
			rev = 'redmax.JointRevolute';
			if ~isa(jointA,rev) || ~isa(jointB,rev)
				error('ConstraintMultQ only works between revolute joints');
			end
		end
		
		%%
		function setFactor(this,factor)
			this.factor = factor;
		end
	end
	
	%%
	methods (Access = protected)
		%%
		function [Gr,Grdot,gr,grdot,grddot] = computeJacEqR_(this,Gr,Grdot,gr,grdot,grddot)
			% Constraint: qB = factor*qA
			% factor*qdotA - qdotB = 0
			% Gr = [factor -1]
			% Grdot = 0
			rows = this.idxER;
			colsA = this.jointA.idxR;
			colsB = this.jointB.idxR;
			this.idxQ = [colsA colsB];
			Gr(rows,colsA) = this.factor;
			Gr(rows,colsB) = -1;
			gr(rows) = this.factor*this.jointA.q(1) - this.jointB.q(1);
		end
	end
end
