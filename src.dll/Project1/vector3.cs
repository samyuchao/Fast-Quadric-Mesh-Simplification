//------------------------------------------------------------------------------
// <auto-generated />
//
// This file was automatically generated by SWIG (http://www.swig.org).
// Version 3.0.12
//
// Do not make changes to this file unless you know what you are doing--modify
// the SWIG interface file instead.
//------------------------------------------------------------------------------


public class vector3 : global::System.IDisposable {
  private global::System.Runtime.InteropServices.HandleRef swigCPtr;
  protected bool swigCMemOwn;

  internal vector3(global::System.IntPtr cPtr, bool cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = new global::System.Runtime.InteropServices.HandleRef(this, cPtr);
  }

  internal static global::System.Runtime.InteropServices.HandleRef getCPtr(vector3 obj) {
    return (obj == null) ? new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero) : obj.swigCPtr;
  }

  ~vector3() {
    Dispose();
  }

  public virtual void Dispose() {
    lock(this) {
      if (swigCPtr.Handle != global::System.IntPtr.Zero) {
        if (swigCMemOwn) {
          swigCMemOwn = false;
          cppdllPINVOKE.delete_vector3(swigCPtr);
        }
        swigCPtr = new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero);
      }
      global::System.GC.SuppressFinalize(this);
    }
  }

  public double x {
    set {
      cppdllPINVOKE.vector3_x_set(swigCPtr, value);
    } 
    get {
      double ret = cppdllPINVOKE.vector3_x_get(swigCPtr);
      return ret;
    } 
  }

  public double y {
    set {
      cppdllPINVOKE.vector3_y_set(swigCPtr, value);
    } 
    get {
      double ret = cppdllPINVOKE.vector3_y_get(swigCPtr);
      return ret;
    } 
  }

  public double z {
    set {
      cppdllPINVOKE.vector3_z_set(swigCPtr, value);
    } 
    get {
      double ret = cppdllPINVOKE.vector3_z_get(swigCPtr);
      return ret;
    } 
  }

  public vector3() : this(cppdllPINVOKE.new_vector3(), true) {
  }

}
