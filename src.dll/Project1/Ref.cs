//------------------------------------------------------------------------------
// <auto-generated />
//
// This file was automatically generated by SWIG (http://www.swig.org).
// Version 3.0.12
//
// Do not make changes to this file unless you know what you are doing--modify
// the SWIG interface file instead.
//------------------------------------------------------------------------------


public class Ref : global::System.IDisposable {
  private global::System.Runtime.InteropServices.HandleRef swigCPtr;
  protected bool swigCMemOwn;

  internal Ref(global::System.IntPtr cPtr, bool cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = new global::System.Runtime.InteropServices.HandleRef(this, cPtr);
  }

  internal static global::System.Runtime.InteropServices.HandleRef getCPtr(Ref obj) {
    return (obj == null) ? new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero) : obj.swigCPtr;
  }

  ~Ref() {
    Dispose();
  }

  public virtual void Dispose() {
    lock(this) {
      if (swigCPtr.Handle != global::System.IntPtr.Zero) {
        if (swigCMemOwn) {
          swigCMemOwn = false;
          cppdllPINVOKE.delete_Ref(swigCPtr);
        }
        swigCPtr = new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero);
      }
      global::System.GC.SuppressFinalize(this);
    }
  }

  public int tid {
    set {
      cppdllPINVOKE.Ref_tid_set(swigCPtr, value);
    } 
    get {
      int ret = cppdllPINVOKE.Ref_tid_get(swigCPtr);
      return ret;
    } 
  }

  public int tvertex {
    set {
      cppdllPINVOKE.Ref_tvertex_set(swigCPtr, value);
    } 
    get {
      int ret = cppdllPINVOKE.Ref_tvertex_get(swigCPtr);
      return ret;
    } 
  }

  public Ref() : this(cppdllPINVOKE.new_Ref(), true) {
  }

}
