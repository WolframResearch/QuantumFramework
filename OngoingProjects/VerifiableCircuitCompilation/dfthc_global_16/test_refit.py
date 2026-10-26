import unittest
import numpy as np
import scipy.linalg as la
import refit as r

class GlobalTests(unittest.TestCase):
    def setUp(self):
        self.rng=np.random.default_rng(137)
        a=self.rng.normal(size=(3,3));self.offset=-(a.T@a+np.eye(3))
        x=self.rng.normal(size=(4,6));self.target=x.T@x
        self.model=r.GlobalModel(self.target,self.offset,7.,(2,3,2))
        self.raw=np.concatenate([self.rng.normal(size=18),self.rng.normal(size=6),[.2,-.3]])

    def differences(self,function,raw,tolerance=2e-6):
        value,gradient=function(raw)
        self.assertTrue(np.isfinite(value));self.assertTrue(np.isfinite(gradient).all())
        for _ in range(8):
            direction=self.rng.normal(size=raw.size);direction/=np.linalg.norm(direction)
            step=1e-5
            finite=(function(raw+step*direction)[0]-function(raw-step*direction)[0])/(2*step)
            self.assertAlmostEqual(finite,float(gradient@direction),delta=tolerance*max(1.,abs(finite)))

    def test_pack_isometry(self):
        a=self.rng.normal(size=(4,3,3));a=(a+a.swapaxes(-1,-2))/2
        self.assertTrue(np.allclose(r.unpack(r.pack(a),3),a))
        self.assertAlmostEqual(float(np.sum(a*a)),float(np.sum(r.pack(a)**2)))

    def test_polar_pullback(self):
        z=self.rng.normal(size=(6,3));g=self.rng.normal(size=(6,3))
        def objective(raw):
            f,s,v=r.polar(raw.reshape(6,3))
            return float(np.sum(f*g)),r.polar_pullback(f,g,s,v).ravel()
        self.differences(objective,z.ravel())

    def test_global_gradient(self):self.differences(self.model.value_gradient,self.raw)

    def test_frame_offset_invariants(self):
        for _ in range(4):
            raw=self.raw+self.rng.normal(size=self.raw.size)*.1
            u,w,v,cache=self.model.materialize(raw)
            x=np.einsum('rbc,rbi,rbj->rcij',w,u,u)
            z=np.einsum('rc,rcij->ij',v,x)
            self.assertTrue(np.allclose(z,self.offset,rtol=0,atol=1e-12))
            self.assertAlmostEqual(float(np.sum(v*v)),7.,places=12)
            self.assertTrue(np.allclose(cache[0].T@cache[0],np.eye(3),rtol=0,atol=1e-13))

    def test_fixed_gradient(self):
        u,w,v,cache=self.model.materialize(self.raw)
        perturbed=u+self.rng.normal(size=u.shape)*.01
        perturbed/=np.linalg.norm(perturbed,axis=-1,keepdims=True)
        fixed=r.FixedModel(self.model,perturbed,-cache[4]**2)
        raw=np.concatenate([w[:,:,1:].ravel(),self.raw[-2:]])
        self.differences(fixed.value_gradient,raw)
        fw,fv,_=fixed.materialize(raw)
        self.assertAlmostEqual(float(np.sum(fv*fv)),7.,places=12)

    def test_global_loss_is_direct(self):
        u,w,_,_=self.model.materialize(self.raw)
        x=np.einsum('rbc,rbd->rcd',w,r.projector_rows(u)).reshape(4,6)
        value,_=self.model.value_gradient(self.raw)
        self.assertAlmostEqual(value*self.model.scale,float(np.sum((x.T@x-self.target)**2)),places=10)

    def test_core_updates(self):
        frame=r.polar(self.raw[:18].reshape(6,3))[0]
        raw,history=r.core_initialize(self.model,frame,3)
        self.assertTrue(np.all(np.diff(history)<=1e-9))
        self.assertAlmostEqual(history[-1],self.model.direct_error(raw),places=10)

    def test_copy_gauge_gradient_and_invariants(self):
        w=self.rng.normal(size=(2,3,3));v=np.zeros((2,3));v[:,0]=[.8,1.1]
        raw=(np.broadcast_to(np.eye(2),(2,2,2))+self.rng.normal(size=(2,2,2))*.1).ravel()
        def objective(q):
            value,gradient,_=r.gauge_value_gradient(q,w,v)
            return value,gradient
        self.differences(objective,raw)
        _,_,rotated=r.gauge_value_gradient(raw,w,v)
        self.assertTrue(np.allclose(w@w.swapaxes(-1,-2),rotated@rotated.swapaxes(-1,-2)))
        self.assertTrue(np.allclose(w@v[:,:,None],rotated@v[:,:,None]))

    def test_domains(self):
        with self.assertRaises(ValueError):r.polar(np.zeros((6,3)))
        with self.assertRaises(ValueError):r.polar(np.ones((2,3)))
        with self.assertRaises(ValueError):r.polar(np.ones((6,3),dtype=complex))
        with self.assertRaises(ValueError):r.GlobalModel(self.target.astype(complex),self.offset,7.,(2,3,2))
        with self.assertRaises(ValueError):r.GlobalModel(self.target,self.offset,np.nan,(2,3,2))
        with self.assertRaises(la.LinAlgError):r.GlobalModel(self.target,np.eye(3),7.,(2,3,2))
        with self.assertRaises(ValueError):self.model.split(np.zeros(5))

if __name__=='__main__':unittest.main()
