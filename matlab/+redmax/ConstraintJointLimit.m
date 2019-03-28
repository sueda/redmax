classdef ConstraintJointLimit < redmax.Constraint
	% ConstraintJointLimit Inequality joint limit constraint
	% Works only for revolute joints
	
	%%
	properties
		joint % Joint to apply the constraint on
		ql    % lower limit
		qu    % upper limit
	end
	
	%%
	methods
		%%
		function this = ConstraintJointLimit(joint)
			this = this@redmax.Constraint(0,0,0,1);
			this.name = [joint.name,'-LIMIT'];
			this.joint = joint;
		end
		
		%%
		function setLimits(this,ql,qu)
			this.ql = ql;
			this.qu = qu;
		end
	end
	
	%%
	methods (Access = protected)
		%%
		function [Cr,Crdot,cr,crdot,crddot] = computeJacIneqR_(this,Cr,Crdot,cr,crdot,crddot)
			rows = this.idxIR;
			cols = this.joint.idxR;
			this.idxQ = cols;
			if this.joint.q(1) <= this.ql
				Cr(rows,cols) = -1;
				cr(rows) = this.ql - this.joint.q(1);
				this.activeR = true;
			elseif this.joint.q(1) >= this.qu
				Cr(rows,cols) = 1;
				cr(rows) = this.qu - this.joint.q(1);
				this.activeR = true;
			else
				this.activeR = false;
			end
		end
		
		%%
		function [value,isterminal,direction] = ineqEventFcn_(this,value,isterminal,direction)
			q = this.joint.q(1);
			value(end+1) = q - this.ql;
			isterminal(end+1) = 1;
			direction(end+1) = -1;
			value(end+1) = q - this.qu;
			isterminal(end+1) = 1;
			direction(end+1) = 1;
		end
		
		%%
		function ineqProjPos_(this)
			if this.joint.q(1) <= this.ql
				this.joint.q(1) = this.ql;
			elseif this.joint.q(1) >= this.qu
				this.joint.q(1) = this.qu;
			end
		end
	end
end
