

#pragma once



#pragma once





#pragma once





#pragma once





#pragma once




#pragma once






#pragma once






 


















































































  



















































































 










































    







    
    


        
            
        


    



















#pragma once





















































































































































































































































































































































































































































































































































































































































































































#pragma region Input Buffer SAL 1 compatibility macros



























































































































































































































































































































































































































































































































































































































































































































































































#pragma endregion Input Buffer SAL 1 compatibility macros
























































































































































































































































































































































































































































































































































































































































































































































































































































































































extern "C" {









































































































































































































































    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    
    
    
    

    
    
























































































































































































































































    
    





















































































}





















#pragma once


extern "C" {











































































































































































































































































































































}


















#pragma once



#pragma pack(push, 8)


extern "C" {







    
    
        typedef unsigned __int64  uintptr_t;
    





    
    


        typedef char* va_list;
    



    













    
    














































    void __cdecl __va_start(va_list* , ...);

    
    



    




} 



    extern "C++"
    {
        template <typename _Ty>
        struct __vcrt_va_list_is_reference
        {
            enum : bool { __the_value = false };
        };

        template <typename _Ty>
        struct __vcrt_va_list_is_reference<_Ty&>
        {
            enum : bool { __the_value = true };
        };

        template <typename _Ty>
        struct __vcrt_va_list_is_reference<_Ty&&>
        {
            enum : bool { __the_value = true };
        };

        template <typename _Ty>
        void __vcrt_va_start_verify_argument_type() throw()
        {
            static_assert(!__vcrt_va_list_is_reference<_Ty>::__the_value, "va_start argument must not have reference type and must not be parenthesized");
        }
    } 

    







#pragma pack(pop)







    



    























__pragma(pack(push, 8)) extern "C" {




    



















    






        
    





    






        
    









    







    
    





    









    







    





    



    
        
        
    









    typedef unsigned __int64 size_t;
    typedef __int64          ptrdiff_t;
    typedef __int64          intptr_t;







    typedef bool  __vcrt_bool;










    



    



    









    
        
    





    





    extern "C++"
    {
        template <typename _CountofType, size_t _SizeOfArray>
        char (*__countof_helper(__unaligned _CountofType (&_Array)[_SizeOfArray]))[_SizeOfArray];

        
    }












    


        




    







    
        
    













    void __cdecl __security_init_cookie(void);

    



        void __cdecl __security_check_cookie(  uintptr_t _StackCookie);
        __declspec(noreturn) void __cdecl __report_gsfailure(  uintptr_t _StackCookie);
    


extern uintptr_t __security_cookie;


    
    
    


} __pragma(pack(pop))











#pragma once



































































































































































































































































































































__pragma(pack(push, 8)) extern "C" {









    


        
    







    



    


        
    







    



    









    






    














    


        
    





    





    





    










extern "C++"
{
    template<bool _Enable, typename _Ty>
    struct _CrtEnableIf;

    template<typename _Ty>
    struct _CrtEnableIf<true, _Ty>
    {
        typedef _Ty _Type;
    };
}



    typedef bool  __crt_bool;













    















    










    











        
    



    



    
        
    




























__declspec(dllimport) void __cdecl _invalid_parameter_noinfo(void);
__declspec(dllimport) __declspec(noreturn) void __cdecl _invalid_parameter_noinfo_noreturn(void);

__declspec(noreturn)
__declspec(dllimport) void __cdecl _invoke_watson(
      wchar_t const*,
      wchar_t const*,
      wchar_t const*,
            unsigned int,
            uintptr_t);


    



        
        
        
        
        
        
        
        
        
        
        
        

    














    


        


    










    






        
    



    


        
    







































    







    





    


        


            
        
    













    


        



    



    
        
    





    
        
        
        
    





    
        
              
        


    





    
        
    





    
        
    







    









typedef int                           errno_t;
typedef unsigned short                wint_t;
typedef unsigned short                wctype_t;
typedef long                          __time32_t;
typedef __int64                       __time64_t;

typedef struct __crt_locale_data_public
{
      unsigned short const* _locale_pctype;
      int _locale_mb_cur_max;
               unsigned int _locale_lc_codepage;
} __crt_locale_data_public;

typedef struct __crt_locale_pointers
{
    struct __crt_locale_data*    locinfo;
    struct __crt_multibyte_data* mbcinfo;
} __crt_locale_pointers;

typedef __crt_locale_pointers* _locale_t; 

typedef struct _Mbstatet
{ 
    unsigned long _Wchar;
    unsigned short _Byte, _State;
} _Mbstatet;

typedef _Mbstatet mbstate_t;










    


        typedef __time64_t time_t;
    




    



    typedef size_t rsize_t;











    

        










        










        










        










        










        










        










        










        










        















        















        
















    




























































































    







































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































        
        
        
        

        


        


        


        


        


        


        


        


        



        



        


        


        


        


        


        


        


        


        


        


        



        



        



        


        



        




        

        




        

        




        

        




        

        




        

        




        

        




        

        




        

    




} __pragma(pack(pop))











    




        
    



    


        
            
        


    



    


        
            
        


    



    




        
            
        


    



#pragma pack(push,8)





 
  


   
  
 

 
  
  
 

















 
  
   


    
   
  






 


 
  


   
  
 


 
  


   
  
 


 
  


   
  
 


 
  


   
  
 


 
  


   
  
 


 
































		

	







		
		


			
		
	


		




 
  
 










































	
	






		


			
		
	

	
	




		


			
		
	

	
	




		
	







	
		#pragma detect_mismatch("_MSC_VER", "1900")
	

	
		#pragma detect_mismatch("_ITERATOR_DEBUG_LEVEL", "0")
	

	
		




			#pragma detect_mismatch("RuntimeLibrary", "MD_DynamicRelease")
		


	









	
		
	
































 


 
 

 









 









 









 







































 
 


 
 

 



























#pragma once























    
    



    










    


#pragma comment(lib, "msvcprt" "" "")






















 















 
  


   
  
 

 
  


   
  
 

 
  


   
  
 


 
  







   



    
   

  
 

 
  
 

 
  


   


     
   
  
 

 


























  
   
  
 

		

 
  
  
  




  
  
  

  







   
   
   
  

  
  
  
  

 














		





		







typedef long long _Longlong;
typedef unsigned long long _ULonglong;

		






		
		






 
namespace std {
enum _Uninitialized
	{	
	_Noinit
	};

		

#pragma warning(push)
#pragma warning(disable:4412)
class __declspec(dllimport) _Lockit
	{	
public:
 

  

















	__thiscall _Lockit();	
	explicit __thiscall _Lockit(int);	
	__thiscall ~_Lockit() noexcept;	
  

	static  void __cdecl _Lockit_ctor(int);
	static  void __cdecl _Lockit_dtor(int);

private:
	static  void __cdecl _Lockit_ctor(_Lockit *);
	static  void __cdecl _Lockit_ctor(_Lockit *, int);
	static  void __cdecl _Lockit_dtor(_Lockit *);

public:
	 _Lockit(const _Lockit&) = delete;
	_Lockit&  operator=(const _Lockit&) = delete;

private:
	int _Locktype;

  












	};

 



































































  



  


  



  


  
 

class __declspec(dllimport) _Init_locks
	{	
public:
 
  











	__thiscall _Init_locks();
	__thiscall ~_Init_locks() noexcept;
  

private:
	static  void __cdecl _Init_locks_ctor(_Init_locks *);
	static  void __cdecl _Init_locks_dtor(_Init_locks *);

 








	};

#pragma warning(pop)
}
 





		

__declspec(dllimport) void __cdecl _Atexit(void (__cdecl *)(void));

typedef unsigned long _Uint32t;




 
 #pragma pack(pop)















 







#pragma once








































































































































































































































































































































typedef signed char        int8_t;
typedef short              int16_t;
typedef int                int32_t;
typedef long long          int64_t;
typedef unsigned char      uint8_t;
typedef unsigned short     uint16_t;
typedef unsigned int       uint32_t;
typedef unsigned long long uint64_t;

typedef signed char        int_least8_t;
typedef short              int_least16_t;
typedef int                int_least32_t;
typedef long long          int_least64_t;
typedef unsigned char      uint_least8_t;
typedef unsigned short     uint_least16_t;
typedef unsigned int       uint_least32_t;
typedef unsigned long long uint_least64_t;

typedef signed char        int_fast8_t;
typedef int                int_fast16_t;
typedef int                int_fast32_t;
typedef long long          int_fast64_t;
typedef unsigned char      uint_fast8_t;
typedef unsigned int       uint_fast16_t;
typedef unsigned int       uint_fast32_t;
typedef unsigned long long uint_fast64_t;

typedef long long          intmax_t;
typedef unsigned long long uintmax_t;










































    
    
    














    



































 
namespace std {
using :: int8_t; using :: int16_t;
using :: int32_t; using :: int64_t;
using :: uint8_t; using :: uint16_t;
using :: uint32_t; using :: uint64_t;

using :: int_least8_t; using :: int_least16_t;
using :: int_least32_t;  using :: int_least64_t;
using :: uint_least8_t; using :: uint_least16_t;
using :: uint_least32_t; using :: uint_least64_t;

using :: int_fast8_t; using :: int_fast16_t;
using :: int_fast32_t;  using :: int_fast64_t;
using :: uint_fast8_t; using :: uint_fast16_t;
using :: uint_fast32_t; using :: uint_fast64_t;

using :: intmax_t; using :: intptr_t;
using :: uintmax_t; using :: uintptr_t;


	namespace tr1 {
using :: int8_t; using :: int16_t;
using :: int32_t; using :: int64_t;
using :: uint8_t; using :: uint16_t;
using :: uint32_t; using :: uint64_t;

using :: int_least8_t; using :: int_least16_t;
using :: int_least32_t;  using :: int_least64_t;
using :: uint_least8_t; using :: uint_least16_t;
using :: uint_least32_t; using :: uint_least64_t;

using :: int_fast8_t; using :: int_fast16_t;
using :: int_fast32_t;  using :: int_fast64_t;
using :: uint_fast8_t; using :: uint_fast16_t;
using :: uint_fast32_t; using :: uint_fast64_t;

using :: intmax_t; using :: intptr_t;
using :: uintmax_t; using :: uintptr_t;
	}	

}
 










#pragma once










 







#pragma once












#pragma once



__pragma(pack(push, 8)) extern "C" {







































     
__declspec(dllimport) __declspec(allocator) __declspec(restrict)
void* __cdecl _calloc_base(
      size_t _Count,
      size_t _Size
    );

     
__declspec(dllimport)  __declspec(allocator) __declspec(restrict)
void* __cdecl calloc(
       size_t _Count,
       size_t _Size
    );

 
__declspec(dllimport) int __cdecl _callnewh(
      size_t _Size
    );

     
__declspec(dllimport) __declspec(allocator)
void* __cdecl _expand(
                void*  _Block,
       size_t _Size
    );

__declspec(dllimport)
void __cdecl _free_base(
        void* _Block
    );

__declspec(dllimport)
void __cdecl free(
        void* _Block
    );

     
__declspec(dllimport) __declspec(allocator) __declspec(restrict)
void* __cdecl _malloc_base(
      size_t _Size
    );

     
__declspec(dllimport) __declspec(allocator)  __declspec(restrict)
void* __cdecl malloc(
       size_t _Size
    );

 
__declspec(dllimport)
size_t __cdecl _msize(
      void* _Block
    );

       
__declspec(dllimport) __declspec(allocator) __declspec(restrict)
void* __cdecl _realloc_base(
         void*  _Block,
                                 size_t _Size
    );

       
__declspec(dllimport) __declspec(allocator) __declspec(restrict)
void* __cdecl realloc(
        void*  _Block,
              size_t _Size
    );

       
__declspec(dllimport) __declspec(allocator) __declspec(restrict)
void* __cdecl _recalloc(
        void*  _Block,
              size_t _Count,
              size_t _Size
    );

__declspec(dllimport)
void __cdecl _aligned_free(
        void* _Block
    );

     
__declspec(dllimport) __declspec(allocator) __declspec(restrict)
void* __cdecl _aligned_malloc(
       size_t _Size,
                         size_t _Alignment
    );

     
__declspec(dllimport) __declspec(allocator) __declspec(restrict)
void* __cdecl _aligned_offset_malloc(
       size_t _Size,
                         size_t _Alignment,
                         size_t _Offset
    );

 
__declspec(dllimport)
size_t __cdecl _aligned_msize(
      void*  _Block,
               size_t _Alignment,
               size_t _Offset
    );

       
__declspec(dllimport) __declspec(allocator) __declspec(restrict)
void* __cdecl _aligned_offset_realloc(
        void*  _Block,
              size_t _Size,
                                size_t _Alignment,
                                size_t _Offset
    );

       
__declspec(dllimport) __declspec(allocator) __declspec(restrict)
void* __cdecl _aligned_offset_recalloc(
        void*  _Block,
              size_t _Count,
              size_t _Size,
                                size_t _Alignment,
                                size_t _Offset
    );

       
__declspec(dllimport) __declspec(allocator) __declspec(restrict)
void* __cdecl _aligned_realloc(
        void*  _Block,
              size_t _Size,
                                size_t _Alignment
    );

       
__declspec(dllimport) __declspec(allocator) __declspec(restrict)
void* __cdecl _aligned_recalloc(
        void*  _Block,
              size_t _Count,
              size_t _Size,
                                size_t _Alignment
    );






















} __pragma(pack(pop))












#pragma once










#pragma once




__pragma(pack(push, 8)) extern "C" {




    namespace std
    {
        typedef decltype(__nullptr) nullptr_t;
    }

    using ::std::nullptr_t;





__declspec(dllimport) int* __cdecl _errno(void);


__declspec(dllimport) errno_t __cdecl _set_errno(  int _Value);
__declspec(dllimport) errno_t __cdecl _get_errno(  int* _Value);



    
        
    






__declspec(dllimport) extern unsigned long  __cdecl __threadid(void);

__declspec(dllimport) extern uintptr_t __cdecl __threadhandle(void);



} __pragma(pack(pop))


__pragma(pack(push, 8)) extern "C" {





     
    __declspec(dllimport) void* __cdecl bsearch_s(
                                                        void const* _Key,
          void const* _Base,
                                                        rsize_t     _NumOfElements,
                                                        rsize_t     _SizeOfElements,
          int (__cdecl* _PtFuncCompare)(void*, void const*, void const*),
                                                        void*       _Context
        );

    __declspec(dllimport) void __cdecl qsort_s(
          void*   _Base,
                                                             rsize_t _NumOfElements,
                                                             rsize_t _SizeOfElements,
          int (__cdecl* _PtFuncCompare)(void*, void const*, void const*),
                                                             void*   _Context
        );





 
__declspec(dllimport) void* __cdecl bsearch(
                                                    void const* _Key,
      void const* _Base,
                                                    size_t      _NumOfElements,
                                                    size_t      _SizeOfElements,
      int (__cdecl* _PtFuncCompare)(void const*, void const*)
    );

__declspec(dllimport) void __cdecl qsort(
      void*  _Base,
                                                         size_t _NumOfElements,
                                                         size_t _SizeOfElements,
      int (__cdecl* _PtFuncCompare)(void const*, void const*)
    );

 
__declspec(dllimport) void* __cdecl _lfind_s(
                                                       void const*   _Key,
      void const*   _Base,
                                                    unsigned int* _NumOfElements,
                                                       size_t        _SizeOfElements,
      int (__cdecl* _PtFuncCompare)(void*, void const*, void const*), 
                                                       void*         _Context
    );

 
__declspec(dllimport) void* __cdecl _lfind(
                                                       void const*   _Key,
      void const*   _Base,
                                                    unsigned int* _NumOfElements,
                                                       unsigned int  _SizeOfElements,
      int (__cdecl* _PtFuncCompare)(void const*, void const*)
    );

 
__declspec(dllimport) void* __cdecl _lsearch_s(
                                                             void const*   _Key,
      void*         _Base,
                                                          unsigned int* _NumOfElements,
                                                             size_t        _SizeOfElements,
      int (__cdecl* _PtFuncCompare)(void*, void const*, void const*),
                                                             void*         _Context
    );

 
__declspec(dllimport) void* __cdecl _lsearch(
                                                             void const*   _Key,
      void*         _Base,
                                                          unsigned int* _NumOfElements,
                                                             unsigned int  _SizeOfElements,
      int (__cdecl* _PtFuncCompare)(void const*, void const*)
    );























































































      __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_lfind" ". See online help for details."))
    __declspec(dllimport) void* __cdecl lfind(
                                                           void const*   _Key,
          void const*   _Base,
                                                        unsigned int* _NumOfElements,
                                                           unsigned int  _SizeOfElements,
          int (__cdecl* _PtFuncCompare)(void const*, void const*)
        );

      __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_lsearch" ". See online help for details."))
    __declspec(dllimport) void* __cdecl lsearch(
                                                                void const*   _Key,
          void*         _Base,
                                                             unsigned int* _NumOfElements,
                                                                unsigned int  _SizeOfElements,
          int (__cdecl* _PtFuncCompare)(void const*, void const*)
        );





} __pragma(pack(pop))










#pragma once



__pragma(pack(push, 8)) extern "C" {
































 

__declspec(dllimport) errno_t __cdecl _itow_s(
                              int      _Value,
      wchar_t* _Buffer,
                              size_t   _BufferCount,
                              int      _Radix
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl _itow_s(  int _Value, wchar_t (&_Buffer)[_Size],   int _Radix) throw() { return _itow_s(_Value, _Buffer, _Size, _Radix); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "_itow_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport) wchar_t* __cdecl _itow( int _Value,   wchar_t *_Buffer,  int _Radix);

 

__declspec(dllimport) errno_t __cdecl _ltow_s(
                              long     _Value,
      wchar_t* _Buffer,
                              size_t   _BufferCount,
                              int      _Radix
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl _ltow_s(  long _Value, wchar_t (&_Buffer)[_Size],   int _Radix) throw() { return _ltow_s(_Value, _Buffer, _Size, _Radix); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "_ltow_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport) wchar_t* __cdecl _ltow( long _Value,   wchar_t *_Buffer,  int _Radix);


__declspec(dllimport) errno_t __cdecl _ultow_s(
                              unsigned long _Value,
      wchar_t*      _Buffer,
                              size_t        _BufferCount,
                              int           _Radix
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl _ultow_s(  unsigned long _Value, wchar_t (&_Buffer)[_Size],   int _Radix) throw() { return _ultow_s(_Value, _Buffer, _Size, _Radix); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "_ultow_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport) wchar_t* __cdecl _ultow( unsigned long _Value,   wchar_t *_Buffer,  int _Radix);

 
__declspec(dllimport) double __cdecl wcstod(
                        wchar_t const* _String,
        wchar_t**      _EndPtr
    );

 
__declspec(dllimport) double __cdecl _wcstod_l(
                        wchar_t const* _String,
        wchar_t**      _EndPtr,
                      _locale_t      _Locale
    );

 
__declspec(dllimport) long __cdecl wcstol(
                        wchar_t const* _String,
        wchar_t**      _EndPtr,
                          int            _Radix
    );

 
__declspec(dllimport) long __cdecl _wcstol_l(
                        wchar_t const* _String,
        wchar_t**      _EndPtr,
                          int            _Radix,
                      _locale_t      _Locale
    );

 
__declspec(dllimport) long long __cdecl wcstoll(
                        wchar_t const* _String,
        wchar_t**      _EndPtr,
                          int            _Radix
    );

 
__declspec(dllimport) long long __cdecl _wcstoll_l(
                        wchar_t const* _String,
        wchar_t**      _EndPtr,
                          int            _Radix,
                      _locale_t      _Locale
    );

 
__declspec(dllimport) unsigned long __cdecl wcstoul(
                        wchar_t const* _String,
        wchar_t**      _EndPtr,
                          int            _Radix
    );

 
__declspec(dllimport) unsigned long __cdecl _wcstoul_l(
                        wchar_t const* _String,
        wchar_t**      _EndPtr,
                          int            _Radix,
                      _locale_t      _Locale
    );

 
__declspec(dllimport) unsigned long long __cdecl wcstoull(
                        wchar_t const* _String,
        wchar_t**      _EndPtr,
                          int            _Radix
    );

 
__declspec(dllimport) unsigned long long __cdecl _wcstoull_l(
                        wchar_t const* _String,
        wchar_t**      _EndPtr,
                          int            _Radix,
                      _locale_t      _Locale
    );

 
__declspec(dllimport) long double __cdecl wcstold(
                        wchar_t const* _String,
        wchar_t**      _EndPtr
    );

 
__declspec(dllimport) long double __cdecl _wcstold_l(
                        wchar_t const* _String,
        wchar_t**      _EndPtr,
                      _locale_t      _Locale
    );

 
__declspec(dllimport) float __cdecl wcstof(
                        wchar_t const* _String,
        wchar_t**      _EndPtr
    );

 
__declspec(dllimport) float __cdecl _wcstof_l(
                        wchar_t const* _String,
        wchar_t**      _EndPtr,
                      _locale_t      _Locale
    );

 
__declspec(dllimport) double __cdecl _wtof(
      wchar_t const* _String
    );

 
__declspec(dllimport) double __cdecl _wtof_l(
        wchar_t const* _String,
      _locale_t      _Locale
    );

 
__declspec(dllimport) int __cdecl _wtoi(
      wchar_t const* _String
    );

 
__declspec(dllimport) int __cdecl _wtoi_l(
        wchar_t const* _String,
      _locale_t      _Locale
    );

 
__declspec(dllimport) long __cdecl _wtol(
      wchar_t const* _String
    );

 
__declspec(dllimport) long __cdecl _wtol_l(
        wchar_t const* _String,
      _locale_t      _Locale
    );

 
__declspec(dllimport) long long __cdecl _wtoll(
      wchar_t const* _String
    );

 
__declspec(dllimport) long long __cdecl _wtoll_l(
        wchar_t const* _String,
      _locale_t      _Locale
    );


__declspec(dllimport) errno_t __cdecl _i64tow_s(
                              __int64  _Value,
      wchar_t* _Buffer,
                              size_t   _BufferCount,
                              int      _Radix
    );

__declspec(deprecated("This function or variable may be unsafe. Consider using " "_i64tow_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__declspec(dllimport) wchar_t* __cdecl _i64tow(
                        __int64  _Value,
        wchar_t* _Buffer,
                        int      _Radix
    );


__declspec(dllimport) errno_t __cdecl _ui64tow_s(
                              unsigned __int64 _Value,
      wchar_t*         _Buffer,
                              size_t           _BufferCount,
                              int              _Radix
    );

__declspec(deprecated("This function or variable may be unsafe. Consider using " "_ui64tow_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__declspec(dllimport) wchar_t* __cdecl _ui64tow(
                        unsigned __int64 _Value,
        wchar_t*         _Buffer,
                        int              _Radix
    );

 
__declspec(dllimport) __int64 __cdecl _wtoi64(
      wchar_t const* _String
    );

 
__declspec(dllimport) __int64 __cdecl _wtoi64_l(
        wchar_t const* _String,
      _locale_t      _Locale
    );

 
__declspec(dllimport) __int64 __cdecl _wcstoi64(
                        wchar_t const* _String,
        wchar_t**      _EndPtr,
                          int            _Radix
    );

 
__declspec(dllimport) __int64 __cdecl _wcstoi64_l(
                        wchar_t const* _String,
        wchar_t**      _EndPtr,
                          int            _Radix,
                      _locale_t      _Locale
    );

 
__declspec(dllimport) unsigned __int64 __cdecl _wcstoui64(
                        wchar_t const* _String,
        wchar_t**      _EndPtr,
                          int            _Radix
    );

 
__declspec(dllimport) unsigned __int64 __cdecl _wcstoui64_l(
                        wchar_t const* _String,
        wchar_t**      _EndPtr,
                          int            _Radix,
                      _locale_t      _Locale
    );




 
 
__declspec(dllimport) __declspec(allocator) wchar_t* __cdecl _wfullpath(
      wchar_t*       _Buffer,
                                wchar_t const* _Path,
                                  size_t         _BufferCount
    );




__declspec(dllimport) errno_t __cdecl _wmakepath_s(
      wchar_t*       _Buffer,
                              size_t         _BufferCount,
                        wchar_t const* _Drive,
                        wchar_t const* _Dir,
                        wchar_t const* _Filename,
                        wchar_t const* _Ext
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl _wmakepath_s(wchar_t (&_Buffer)[_Size],   wchar_t const* _Drive,   wchar_t const* _Dir,   wchar_t const* _Filename,   wchar_t const* _Ext) throw() { return _wmakepath_s(_Buffer, _Size, _Drive, _Dir, _Filename, _Ext); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "_wmakepath_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport) void __cdecl _wmakepath(  wchar_t *_Buffer,  wchar_t const* _Drive,  wchar_t const* _Dir,  wchar_t const* _Filename,  wchar_t const* _Ext);

__declspec(dllimport) void __cdecl _wperror(
      wchar_t const* _ErrMsg
    );

__declspec(deprecated("This function or variable may be unsafe. Consider using " "_wsplitpath_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__declspec(dllimport) void __cdecl _wsplitpath(
                        wchar_t const* _FullPath,
        wchar_t*       _Drive,
        wchar_t*       _Dir,
        wchar_t*       _Filename,
        wchar_t*       _Ext
    );

__declspec(dllimport) errno_t __cdecl _wsplitpath_s(
                                  wchar_t const* _FullPath,
         wchar_t*       _Drive,
                                    size_t         _DriveCount,
           wchar_t*       _Dir,
                                    size_t         _DirCount,
      wchar_t*       _Filename,
                                    size_t         _FilenameCount,
           wchar_t*       _Ext,
                                    size_t         _ExtCount
    );

extern "C++" { template <size_t _DriveSize, size_t _DirSize, size_t _NameSize, size_t _ExtSize> inline errno_t __cdecl _wsplitpath_s(   wchar_t const* _Path,   wchar_t (&_Drive)[_DriveSize],   wchar_t (&_Dir)[_DirSize],   wchar_t (&_Name)[_NameSize],   wchar_t (&_Ext)[_ExtSize] ) throw() { return _wsplitpath_s(_Path, _Drive, _DriveSize, _Dir, _DirSize, _Name, _NameSize, _Ext, _ExtSize); } }





    
    

    
    __declspec(dllimport) errno_t __cdecl _wdupenv_s(
            wchar_t**      _Buffer,
                                                                            size_t*        _BufferCount,
                                                                               wchar_t const* _VarName
        );

    

      __declspec(deprecated("This function or variable may be unsafe. Consider using " "_wdupenv_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
    __declspec(dllimport) wchar_t* __cdecl _wgetenv(
          wchar_t const* _VarName
        );

     
    
    __declspec(dllimport) errno_t __cdecl _wgetenv_s(
                                     size_t*        _RequiredCount,
          wchar_t*       _Buffer,
                                      size_t         _BufferCount,
                                    wchar_t const* _VarName
        );

    extern "C++" { template <size_t _Size> inline   errno_t __cdecl _wgetenv_s(  size_t* _RequiredCount, wchar_t (&_Buffer)[_Size],   wchar_t const* _VarName) throw() { return _wgetenv_s(_RequiredCount, _Buffer, _Size, _VarName); } }

     
    __declspec(dllimport) int __cdecl _wputenv(
          wchar_t const* _EnvString
        );

    
    __declspec(dllimport) errno_t __cdecl _wputenv_s(
          wchar_t const* _Name,
          wchar_t const* _Value
        );
    
    __declspec(dllimport) errno_t __cdecl _wsearchenv_s(
                                wchar_t const* _Filename,
                                wchar_t const* _VarName,
          wchar_t*       _Buffer,
                                  size_t         _BufferCount
        );
    
    extern "C++" { template <size_t _Size> inline errno_t __cdecl _wsearchenv_s(  wchar_t const* _Filename,   wchar_t const* _VarName, wchar_t (&_ResultPath)[_Size]) throw() { return _wsearchenv_s(_Filename, _VarName, _ResultPath, _Size); } }
    
    __declspec(deprecated("This function or variable may be unsafe. Consider using " "_wsearchenv_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport) void __cdecl _wsearchenv( wchar_t const* _Filename,  wchar_t const* _VarName,   wchar_t *_ResultPath);

    __declspec(dllimport) int __cdecl _wsystem(
          wchar_t const* _Command
        );





} __pragma(pack(pop))









#pragma once




































































































































































































































































































































__pragma(pack(push, 8)) extern "C" {









    
    
























































    
        
    




} __pragma(pack(pop))


__pragma(pack(push, 8)) extern "C" {




    










__declspec(dllimport) void __cdecl _swab(
        char* _Buf1,
        char* _Buf2,
                                                                  int   _SizeInBytes
    );












__declspec(dllimport) __declspec(noreturn) void __cdecl exit(  int _Code);
__declspec(dllimport) __declspec(noreturn) void __cdecl _exit(  int _Code);
__declspec(dllimport) __declspec(noreturn) void __cdecl _Exit(  int _Code);
__declspec(dllimport) __declspec(noreturn) void __cdecl quick_exit(  int _Code);
__declspec(dllimport) __declspec(noreturn) void __cdecl abort(void);





__declspec(dllimport) unsigned int __cdecl _set_abort_behavior(
      unsigned int _Flags,
      unsigned int _Mask
    );




    typedef int (__cdecl* _onexit_t)(void);










    
    


























































    int       __cdecl atexit(void (__cdecl*)(void));
    _onexit_t __cdecl _onexit(  _onexit_t _Func);


int __cdecl at_quick_exit(void (__cdecl*)(void));









    
    typedef void (__cdecl* _purecall_handler)(void);

    
    typedef void (__cdecl* _invalid_parameter_handler)(
        wchar_t const*,
        wchar_t const*,
        wchar_t const*, 
        unsigned int,
        uintptr_t
        );

    
    __declspec(dllimport) _purecall_handler __cdecl _set_purecall_handler(
          _purecall_handler _Handler
        );

    __declspec(dllimport) _purecall_handler __cdecl _get_purecall_handler(void);

    
    __declspec(dllimport) _invalid_parameter_handler __cdecl _set_invalid_parameter_handler(
          _invalid_parameter_handler _Handler
        );

    __declspec(dllimport) _invalid_parameter_handler __cdecl _get_invalid_parameter_handler(void);

    __declspec(dllimport) _invalid_parameter_handler __cdecl _set_thread_local_invalid_parameter_handler(
          _invalid_parameter_handler _Handler
        );

    __declspec(dllimport) _invalid_parameter_handler __cdecl _get_thread_local_invalid_parameter_handler(void);























 __declspec(dllimport) int __cdecl _set_error_mode(  int _Mode);



__declspec(dllimport) int* __cdecl _errno(void);


__declspec(dllimport) errno_t __cdecl _set_errno(  int _Value);
__declspec(dllimport) errno_t __cdecl _get_errno(  int* _Value);

__declspec(dllimport) unsigned long* __cdecl __doserrno(void);


__declspec(dllimport) errno_t __cdecl _set_doserrno(  unsigned long _Value);
__declspec(dllimport) errno_t __cdecl _get_doserrno(  unsigned long * _Value);


__declspec(dllimport) __declspec(deprecated("This function or variable may be unsafe. Consider using " "strerror" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) char** __cdecl __sys_errlist(void);


__declspec(dllimport) __declspec(deprecated("This function or variable may be unsafe. Consider using " "strerror" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) int * __cdecl __sys_nerr(void);


__declspec(dllimport) void __cdecl perror(  char const* _ErrMsg);




__declspec(deprecated("This function or variable may be unsafe. Consider using " "_get_pgmptr" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport) char**    __cdecl __p__pgmptr (void);
__declspec(deprecated("This function or variable may be unsafe. Consider using " "_get_wpgmptr" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport) wchar_t** __cdecl __p__wpgmptr(void);
__declspec(deprecated("This function or variable may be unsafe. Consider using " "_get_fmode" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport) int*      __cdecl __p__fmode  (void);






    
    
    


 
__declspec(dllimport) errno_t __cdecl _get_pgmptr (  char**    _Value);

 
__declspec(dllimport) errno_t __cdecl _get_wpgmptr(  wchar_t** _Value);

__declspec(dllimport) errno_t __cdecl _set_fmode  (               int       _Mode );

__declspec(dllimport) errno_t __cdecl _get_fmode  (              int*      _PMode);








typedef struct _div_t
{
    int quot;
    int rem;
} div_t;

typedef struct _ldiv_t
{
    long quot;
    long rem;
} ldiv_t;

typedef struct _lldiv_t
{
    long long quot;
    long long rem;
} lldiv_t;

  int       __cdecl abs   (  int       _Number);
  long      __cdecl labs  (  long      _Number);
  long long __cdecl llabs (  long long _Number);
  __int64   __cdecl _abs64(  __int64   _Number);

  unsigned short   __cdecl _byteswap_ushort(  unsigned short   _Number);
  unsigned long    __cdecl _byteswap_ulong (  unsigned long    _Number);
  unsigned __int64 __cdecl _byteswap_uint64(  unsigned __int64 _Number);

  __declspec(dllimport) div_t   __cdecl div  (  int       _Numerator,   int       _Denominator);
  __declspec(dllimport) ldiv_t  __cdecl ldiv (  long      _Numerator,   long      _Denominator);
  __declspec(dllimport) lldiv_t __cdecl lldiv(  long long _Numerator,   long long _Denominator);



#pragma warning (push)
#pragma warning (disable:6540) 

unsigned int __cdecl _rotl(
      unsigned int _Value,
      int          _Shift
    );

 
unsigned long __cdecl _lrotl(
      unsigned long _Value,
      int           _Shift
    );

unsigned __int64 __cdecl _rotl64(
      unsigned __int64 _Value,
      int              _Shift
    );

unsigned int __cdecl _rotr(
      unsigned int _Value,
      int          _Shift
    );

 
unsigned long __cdecl _lrotr(
      unsigned long _Value,
      int           _Shift
    );

unsigned __int64 __cdecl _rotr64(
      unsigned __int64 _Value,
      int              _Shift
    );

#pragma warning (pop)






__declspec(dllimport) void __cdecl srand(  unsigned int _Seed);

  __declspec(dllimport) int __cdecl rand(void);








extern "C++"
{
    inline long abs(long const _X) throw()
    {
        return labs(_X);
    }

    inline long long abs(long long const _X) throw()
    {
        return llabs(_X);
    }

    inline ldiv_t div(long const _A1, long const _A2) throw()
    {
        return ldiv(_A1, _A2);
    }

    inline lldiv_t div(long long const _A1, long long const _A2) throw()
    {
        return lldiv(_A1, _A2);
    }
}











    #pragma pack(push, 4)
    typedef struct
    {
        unsigned char ld[10];
    } _LDOUBLE;
    #pragma pack(pop)

    













typedef struct
{
    double x;
} _CRT_DOUBLE;

typedef struct
{
    float f;
} _CRT_FLOAT;





typedef struct
{
    long double x;
} _LONGDOUBLE;



#pragma pack(push, 4)
typedef struct
{
    unsigned char ld12[12];
} _LDBL12;
#pragma pack(pop)








                     __declspec(dllimport) double    __cdecl atof   (  char const* _String);
   __declspec(dllimport) int       __cdecl atoi   (  char const* _String);
                     __declspec(dllimport) long      __cdecl atol   (  char const* _String);
                     __declspec(dllimport) long long __cdecl atoll  (  char const* _String);
                     __declspec(dllimport) __int64   __cdecl _atoi64(  char const* _String);

  __declspec(dllimport) double    __cdecl _atof_l  (  char const* _String,   _locale_t _Locale);
  __declspec(dllimport) int       __cdecl _atoi_l  (  char const* _String,   _locale_t _Locale);
  __declspec(dllimport) long      __cdecl _atol_l  (  char const* _String,   _locale_t _Locale);
  __declspec(dllimport) long long __cdecl _atoll_l (  char const* _String,   _locale_t _Locale);
  __declspec(dllimport) __int64   __cdecl _atoi64_l(  char const* _String,   _locale_t _Locale);

  __declspec(dllimport) int __cdecl _atoflt (  _CRT_FLOAT*  _Result,   char const* _String);
  __declspec(dllimport) int __cdecl _atodbl (  _CRT_DOUBLE* _Result,   char*       _String);
  __declspec(dllimport) int __cdecl _atoldbl(  _LDOUBLE*    _Result,   char*       _String);

 
__declspec(dllimport) int __cdecl _atoflt_l(
         _CRT_FLOAT* _Result,
        char const* _String,
      _locale_t   _Locale
    );

 
__declspec(dllimport) int __cdecl _atodbl_l(
         _CRT_DOUBLE* _Result,
        char*        _String,
      _locale_t    _Locale
    );


 
__declspec(dllimport) int __cdecl _atoldbl_l(
         _LDOUBLE* _Result,
        char*     _String,
      _locale_t _Locale
    );

 
__declspec(dllimport) float __cdecl strtof(
                        char const* _String,
        char**      _EndPtr
    );

 
__declspec(dllimport) float __cdecl _strtof_l(
                        char const* _String,
        char**      _EndPtr,
                      _locale_t   _Locale
    );

 
__declspec(dllimport) double __cdecl strtod(
                        char const* _String,
        char**      _EndPtr
    );

 
__declspec(dllimport) double __cdecl _strtod_l(
                        char const* _String,
        char**      _EndPtr,
                      _locale_t   _Locale
    );

 
__declspec(dllimport) long double __cdecl strtold(
                        char const* _String,
        char**      _EndPtr
    );

 
__declspec(dllimport) long double __cdecl _strtold_l(
                        char const* _String,
        char**      _EndPtr,
                      _locale_t   _Locale
    );

 
__declspec(dllimport) long __cdecl strtol(
                        char const* _String,
        char**      _EndPtr,
                          int         _Radix
    );

 
__declspec(dllimport) long __cdecl _strtol_l(
                        char const* _String,
        char**      _EndPtr,
                          int         _Radix,
                      _locale_t   _Locale
    );

 
__declspec(dllimport) long long __cdecl strtoll(
                        char const* _String,
        char**      _EndPtr,
                          int         _Radix
    );

 
__declspec(dllimport) long long __cdecl _strtoll_l(
                        char const* _String,
        char**      _EndPtr,
                          int         _Radix,
                      _locale_t   _Locale
    );

 
__declspec(dllimport) unsigned long __cdecl strtoul(
                        char const* _String,
        char**      _EndPtr,
                          int         _Radix
    );

 
__declspec(dllimport) unsigned long __cdecl _strtoul_l(
                        char const* _String,
        char**      _EndPtr,
                          int         _Radix,
                      _locale_t   _Locale
    );

 
__declspec(dllimport) unsigned long long __cdecl strtoull(
                        char const* _String,
        char**      _EndPtr,
                          int         _Radix
    );

 
__declspec(dllimport) unsigned long long __cdecl _strtoull_l(
                        char const* _String,
        char**      _EndPtr,
                          int         _Radix,
                      _locale_t   _Locale
    );

 
__declspec(dllimport) __int64 __cdecl _strtoi64(
                        char const* _String,
        char**      _EndPtr,
                          int         _Radix
    );

 
__declspec(dllimport) __int64 __cdecl _strtoi64_l(
                        char const* _String,
        char**      _EndPtr,
                          int         _Radix,
                      _locale_t   _Locale
    );

 
__declspec(dllimport) unsigned __int64 __cdecl _strtoui64(
                        char const* _String,
        char**      _EndPtr,
                          int         _Radix
    );

 
__declspec(dllimport) unsigned __int64 __cdecl _strtoui64_l(
                        char const* _String,
        char**      _EndPtr,
                          int         _Radix,
                      _locale_t   _Locale
    );








 

__declspec(dllimport) errno_t __cdecl _itoa_s(
                              int    _Value,
      char*  _Buffer,
                              size_t _BufferCount,
                              int    _Radix
    );

extern "C++" { template <size_t _Size> inline   errno_t __cdecl _itoa_s(  int _Value, char (&_Buffer)[_Size],   int _Radix) throw() { return _itoa_s(_Value, _Buffer, _Size, _Radix); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "_itoa_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport) char* __cdecl _itoa( int _Value,   char *_Buffer,  int _Radix);

 

__declspec(dllimport) errno_t __cdecl _ltoa_s(
                              long   _Value,
      char*  _Buffer,
                              size_t _BufferCount,
                              int    _Radix
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl _ltoa_s(  long _Value, char (&_Buffer)[_Size],   int _Radix) throw() { return _ltoa_s(_Value, _Buffer, _Size, _Radix); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "_ltoa_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport) char* __cdecl _ltoa( long _Value,   char *_Buffer,  int _Radix);

 

__declspec(dllimport) errno_t __cdecl _ultoa_s(
                              unsigned long _Value,
      char*         _Buffer,
                              size_t        _BufferCount,
                              int           _Radix
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl _ultoa_s(  unsigned long _Value, char (&_Buffer)[_Size],   int _Radix) throw() { return _ultoa_s(_Value, _Buffer, _Size, _Radix); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "_ultoa_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport) char* __cdecl _ultoa( unsigned long _Value,   char *_Buffer,  int _Radix);

 

__declspec(dllimport) errno_t __cdecl _i64toa_s(
                              __int64 _Value,
      char*   _Buffer,
                              size_t  _BufferCount,
                              int     _Radix
    );

 
__declspec(deprecated("This function or variable may be unsafe. Consider using " "_i64toa_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__declspec(dllimport) char* __cdecl _i64toa(
                        __int64 _Value,
        char*   _Buffer,
                        int     _Radix
    );

 

__declspec(dllimport) errno_t __cdecl _ui64toa_s(
                              unsigned __int64 _Value,
      char*            _Buffer,
                              size_t           _BufferCount,
                              int              _Radix
    );

__declspec(deprecated("This function or variable may be unsafe. Consider using " "_ui64toa_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__declspec(dllimport) char* __cdecl _ui64toa(
                        unsigned __int64 _Value,
        char*            _Buffer,
                        int              _Radix
    );













 

__declspec(dllimport) errno_t __cdecl _ecvt_s(
      char* _Buffer,
       size_t                       _BufferCount,
       double                       _Value,
       int                          _DigitCount,
      int*                         _PtDec,
      int*                         _PtSign
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl _ecvt_s(char (&_Buffer)[_Size],   double _Value,   int _DigitCount,   int* _PtDec,   int* _PtSign) throw() { return _ecvt_s(_Buffer, _Size, _Value, _DigitCount, _PtDec, _PtSign); } }

  __declspec(deprecated("This function or variable may be unsafe. Consider using " "_ecvt_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__declspec(dllimport) char* __cdecl _ecvt(
       double _Value,
       int    _DigitCount,
      int*   _PtDec,
      int*   _PtSign
    );

 

__declspec(dllimport) errno_t __cdecl _fcvt_s(
      char*  _Buffer,
                              size_t _BufferCount,
                              double _Value,
                              int    _FractionalDigitCount,
                             int*   _PtDec,
                             int*   _PtSign
    );

extern "C++" { template <size_t _Size> inline   errno_t __cdecl _fcvt_s(char (&_Buffer)[_Size],   double _Value,   int _FractionalDigitCount,   int* _PtDec,   int* _PtSign) throw() { return _fcvt_s(_Buffer, _Size, _Value, _FractionalDigitCount, _PtDec, _PtSign); } }

 
  __declspec(deprecated("This function or variable may be unsafe. Consider using " "_fcvt_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__declspec(dllimport) char* __cdecl _fcvt(
       double _Value,
       int    _FractionalDigitCount,
      int*   _PtDec,
      int*   _PtSign
    );

 
__declspec(dllimport) errno_t __cdecl _gcvt_s(
      char*  _Buffer,
                              size_t _BufferCount,
                              double _Value,
                              int    _DigitCount
    );

extern "C++" { template <size_t _Size> inline   errno_t __cdecl _gcvt_s(char (&_Buffer)[_Size],   double _Value,   int _DigitCount) throw() { return _gcvt_s(_Buffer, _Size, _Value, _DigitCount); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "_gcvt_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__declspec(dllimport) char* __cdecl _gcvt(
                        double _Value,
                        int    _DigitCount,
        char*  _Buffer
    );











    


        
    

    


        
    

     
    __declspec(dllimport) int __cdecl ___mb_cur_max_func(void);

     
    __declspec(dllimport) int __cdecl ___mb_cur_max_l_func(_locale_t);




 
__declspec(dllimport) int __cdecl mblen(
        char const* _Ch,
                                             size_t      _MaxCount
    );

  
__declspec(dllimport) int __cdecl _mblen_l(
        char const* _Ch,
                                             size_t      _MaxCount,
                                         _locale_t   _Locale
    );

 
 
__declspec(dllimport) size_t __cdecl _mbstrlen(
      char const* _String
    );

 
 
__declspec(dllimport) size_t __cdecl _mbstrlen_l(
        char const* _String, 
      _locale_t   _Locale
    );

 
 
__declspec(dllimport) size_t __cdecl _mbstrnlen(
      char const* _String,
        size_t      _MaxCount
    );

 
 
__declspec(dllimport) size_t __cdecl _mbstrnlen_l(
        char const* _String,
          size_t      _MaxCount,
      _locale_t   _Locale
    );

 
__declspec(dllimport) int __cdecl mbtowc(
                      wchar_t*    _DstCh,
      char const* _SrcCh,
                                      size_t      _SrcSizeInBytes
    );

 
__declspec(dllimport) int __cdecl _mbtowc_l(
                      wchar_t*    _DstCh,
      char const* _SrcCh,
                                      size_t      _SrcSizeInBytes,
                                  _locale_t   _Locale
    );


__declspec(dllimport) errno_t __cdecl mbstowcs_s(
                                                      size_t*     _PtNumOfCharConverted,
      wchar_t*    _DstBuf,
                                                           size_t      _SizeInWords,
                                     char const* _SrcBuf,
                                                           size_t      _MaxCount
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl mbstowcs_s(  size_t* _PtNumOfCharConverted,   wchar_t (&_Dest)[_Size],   char const* _Source,   size_t _MaxCount) throw() { return mbstowcs_s(_PtNumOfCharConverted, _Dest, _Size, _Source, _MaxCount); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "mbstowcs_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport) size_t __cdecl mbstowcs( wchar_t *_Dest,  char const* _Source,  size_t _MaxCount);


__declspec(dllimport) errno_t __cdecl _mbstowcs_s_l(
                                                      size_t*     _PtNumOfCharConverted,
      wchar_t*    _DstBuf,
                                                           size_t      _SizeInWords,
                                     char const* _SrcBuf,
                                                           size_t      _MaxCount,
                                                       _locale_t   _Locale
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl _mbstowcs_s_l(  size_t* _PtNumOfCharConverted,   wchar_t (&_Dest)[_Size],   char const* _Source,   size_t _MaxCount,   _locale_t _Locale) throw() { return _mbstowcs_s_l(_PtNumOfCharConverted, _Dest, _Size, _Source, _MaxCount, _Locale); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "_mbstowcs_s_l" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport) size_t __cdecl _mbstowcs_l(  wchar_t *_Dest,   char const* _Source,   size_t _MaxCount,   _locale_t _Locale);




__declspec(deprecated("This function or variable may be unsafe. Consider using " "wctomb_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__declspec(dllimport) int __cdecl wctomb(
      char*   _MbCh,
                                wchar_t _WCh
    );

__declspec(deprecated("This function or variable may be unsafe. Consider using " "_wctomb_s_l" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__declspec(dllimport) int __cdecl _wctomb_l(
        char*     _MbCh,
                          wchar_t   _WCh,
                      _locale_t _Locale
    );



    
    __declspec(dllimport) errno_t __cdecl wctomb_s(
                                                         int*    _SizeConverted,
          char*   _MbCh,
                                                              rsize_t _SizeInBytes,
                                                              wchar_t _WCh
        );




__declspec(dllimport) errno_t __cdecl _wctomb_s_l(
                             int*     _SizeConverted,
      char*     _MbCh,
                                  size_t    _SizeInBytes,
                                  wchar_t   _WCh, 
                              _locale_t _Locale);


__declspec(dllimport) errno_t __cdecl wcstombs_s(
                                                               size_t*        _PtNumOfCharConverted,
      char*          _Dst,
                                                                    size_t         _DstSizeInBytes,
                                                                  wchar_t const* _Src,
                                                                    size_t         _MaxCountInBytes
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl wcstombs_s(  size_t* _PtNumOfCharConverted,   char (&_Dest)[_Size],   wchar_t const* _Source,   size_t _MaxCount) throw() { return wcstombs_s(_PtNumOfCharConverted, _Dest, _Size, _Source, _MaxCount); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "wcstombs_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport) size_t __cdecl wcstombs( char *_Dest,  wchar_t const* _Source,  size_t _MaxCount);


__declspec(dllimport) errno_t __cdecl _wcstombs_s_l(
                                                               size_t*        _PtNumOfCharConverted,
      char*          _Dst,
                                                                    size_t         _DstSizeInBytes,
                                                                  wchar_t const* _Src,
                                                                    size_t         _MaxCountInBytes,
                                                                _locale_t      _Locale
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl _wcstombs_s_l(  size_t* _PtNumOfCharConverted,   char (&_Dest)[_Size],   wchar_t const* _Source,   size_t _MaxCount,   _locale_t _Locale) throw() { return _wcstombs_s_l(_PtNumOfCharConverted, _Dest, _Size, _Source, _MaxCount, _Locale); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "_wcstombs_s_l" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport) size_t __cdecl _wcstombs_l(  char *_Dest,   wchar_t const* _Source,   size_t _MaxCount,   _locale_t _Locale);




















 
 
__declspec(dllimport) __declspec(allocator) char* __cdecl _fullpath(
      char*       _Buffer,
                                char const* _Path,
                                  size_t      _BufferCount
    );




__declspec(dllimport) errno_t __cdecl _makepath_s(
      char*       _Buffer,
                              size_t      _BufferCount,
                        char const* _Drive,
                        char const* _Dir,
                        char const* _Filename,
                        char const* _Ext
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl _makepath_s(char (&_Buffer)[_Size],   char const* _Drive,   char const* _Dir,   char const* _Filename,   char const* _Ext) throw() { return _makepath_s(_Buffer, _Size, _Drive, _Dir, _Filename, _Ext); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "_makepath_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport) void __cdecl _makepath(  char *_Buffer,  char const* _Drive,  char const* _Dir,  char const* _Filename,  char const* _Ext);

__declspec(deprecated("This function or variable may be unsafe. Consider using " "_splitpath_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__declspec(dllimport) void __cdecl _splitpath(
                        char const* _FullPath,
        char*       _Drive,
        char*       _Dir,
        char*       _Filename,
        char*       _Ext
    );


__declspec(dllimport) errno_t __cdecl _splitpath_s(
                                  char const* _FullPath,
         char*       _Drive,
                                    size_t      _DriveCount,
           char*       _Dir,
                                    size_t      _DirCount,
      char*       _Filename,
                                    size_t      _FilenameCount,
           char*       _Ext,
                                    size_t      _ExtCount
    );

extern "C++" { template <size_t _DriveSize, size_t _DirSize, size_t _NameSize, size_t _ExtSize> inline errno_t __cdecl _splitpath_s(   char const* _Dest,   char (&_Drive)[_DriveSize],   char (&_Dir)[_DirSize],   char (&_Name)[_NameSize],   char (&_Ext)[_ExtSize] ) throw() { return _splitpath_s(_Dest, _Drive, _DriveSize, _Dir, _DirSize, _Name, _NameSize, _Ext, _ExtSize); } }










    

    
     
    __declspec(dllimport) errno_t __cdecl getenv_s(
                                     size_t*     _RequiredCount,
          char*       _Buffer,
                                      rsize_t     _BufferCount,
                                    char const* _VarName
        );
    
    




    __declspec(dllimport) int*       __cdecl __p___argc (void);
    __declspec(dllimport) char***    __cdecl __p___argv (void);
    __declspec(dllimport) wchar_t*** __cdecl __p___wargv(void);

    




        
        
        
    
    
    __declspec(dllimport) char***    __cdecl __p__environ (void);
    __declspec(dllimport) wchar_t*** __cdecl __p__wenviron(void);

    
        
    
    
    





        
        
    



    
    



      __declspec(deprecated("This function or variable may be unsafe. Consider using " "_dupenv_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) 
    __declspec(dllimport) char* __cdecl getenv(
          char const* _VarName
        );

    extern "C++" { template <size_t _Size> inline errno_t __cdecl getenv_s(  size_t* _RequiredCount, char (&_Buffer)[_Size],   char const* _VarName) throw() { return getenv_s(_RequiredCount, _Buffer, _Size, _VarName); } }

    




    
    __declspec(dllimport) errno_t __cdecl _dupenv_s(
            char**      _Buffer,
                                                                            size_t*     _BufferCount,
                                                                               char const* _VarName
        );

    



    __declspec(dllimport) int __cdecl system(
          char const* _Command
        );

    
    
    #pragma warning (push)
    #pragma warning (disable:6540)

     
    __declspec(dllimport) int __cdecl _putenv(
          char const* _EnvString
        );

    
    __declspec(dllimport) errno_t __cdecl _putenv_s(
          char const* _Name,
          char const* _Value
        );

    #pragma warning (pop)

    __declspec(dllimport) errno_t __cdecl _searchenv_s(
                                char const* _Filename,
                                char const* _VarName,
          char*       _Buffer,
                                  size_t      _BufferCount
        );

    extern "C++" { template <size_t _Size> inline errno_t __cdecl _searchenv_s(  char const* _Filename,   char const* _VarName, char (&_Buffer)[_Size]) throw() { return _searchenv_s(_Filename, _VarName, _Buffer, _Size); } }

    __declspec(deprecated("This function or variable may be unsafe. Consider using " "_searchenv_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport) void __cdecl _searchenv( char const* _Filename,  char const* _VarName,   char *_Buffer);

    
    __declspec(deprecated("This function or variable has been superceded by newer library " "or operating system functionality. Consider using " "SetErrorMode" " " "instead. See online help for details."))
    __declspec(dllimport) void __cdecl _seterrormode(
          int _Mode
        );

    __declspec(deprecated("This function or variable has been superceded by newer library " "or operating system functionality. Consider using " "Beep" " " "instead. See online help for details."))
    __declspec(dllimport) void __cdecl _beep(
          unsigned _Frequency,
          unsigned _Duration
        );

    __declspec(deprecated("This function or variable has been superceded by newer library " "or operating system functionality. Consider using " "Sleep" " " "instead. See online help for details."))
    __declspec(dllimport) void __cdecl _sleep(
          unsigned long _Duration
        );












    




    
    

    #pragma warning(push)
    #pragma warning(disable: 4141) 

      __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_ecvt" ". See online help for details.")) __declspec(deprecated("This function or variable may be unsafe. Consider using " "_ecvt_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
    __declspec(dllimport) char* __cdecl ecvt(
           double _Value,
           int    _DigitCount,
          int*   _PtDec,
          int*   _PtSign
        );

      __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_fcvt" ". See online help for details.")) __declspec(deprecated("This function or variable may be unsafe. Consider using " "_fcvt_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
    __declspec(dllimport) char* __cdecl fcvt(
           double _Value,
           int    _FractionalDigitCount,
          int*   _PtDec,
          int*   _PtSign
        );

    __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_gcvt" ". See online help for details.")) __declspec(deprecated("This function or variable may be unsafe. Consider using " "_fcvt_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
    __declspec(dllimport) char* __cdecl gcvt(
                            double _Value,
                            int    _DigitCount,
            char*  _DstBuf
        );
    
    __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_itoa" ". See online help for details.")) __declspec(deprecated("This function or variable may be unsafe. Consider using " "_itoa_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
    __declspec(dllimport) char* __cdecl itoa(
                            int   _Value,
            char* _Buffer,
                            int   _Radix
        );
    
    __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_ltoa" ". See online help for details.")) __declspec(deprecated("This function or variable may be unsafe. Consider using " "_ltoa_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
    __declspec(dllimport) char* __cdecl ltoa(
                            long  _Value,
            char* _Buffer,
                            int   _Radix
        );


    __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_swab" ". See online help for details."))
    __declspec(dllimport) void __cdecl swab(
          char* _Buf1,
          char* _Buf2,
                                     int   _SizeInBytes
        );

    __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_ultoa" ". See online help for details.")) __declspec(deprecated("This function or variable may be unsafe. Consider using " "_ultoa_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
    __declspec(dllimport) char* __cdecl ultoa(
                            unsigned long _Value,
            char*         _Buffer,
                            int           _Radix
        );

    

        

          __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_putenv" ". See online help for details."))
        __declspec(dllimport) int __cdecl putenv(
              char const* _EnvString
            );

    

    #pragma warning(pop)

    _onexit_t __cdecl onexit(  _onexit_t _Func);





} __pragma(pack(pop))




 
namespace std {
using :: size_t; using :: div_t; using :: ldiv_t;
using :: abort; using :: abs; using :: atexit;
using :: atof; using :: atoi; using :: atol;
using :: bsearch; using :: calloc; using :: div;
using :: exit; using :: free;
using :: labs; using :: ldiv; using :: malloc;
using :: mblen; using :: mbstowcs; using :: mbtowc;
using :: qsort; using :: rand; using :: realloc;
using :: srand; using :: strtod; using :: strtol;
using :: strtoul;
using :: wcstombs; using :: wctomb;

using :: lldiv_t;

 
using :: getenv;
using :: system;
 

using :: atoll; using :: llabs; using :: lldiv;
using :: strtof; using :: strtold;
using :: strtoll; using :: strtoull;

using :: _Exit; using :: at_quick_exit; using :: quick_exit;
}
 










#pragma once





#pragma once





 #pragma pack(push,8)
 #pragma warning(push,3)
 

  

 
 
extern "C" {
 
 

		





		






void __cdecl _Feraise(int);

typedef union
	{	
	unsigned short _Word[8];
	float _Float;
	double _Double;
	long double _Long_double;
	} _Dconst;

		
__declspec(dllimport) double __cdecl _Cosh(double, double);
__declspec(dllimport) short __cdecl _Dtest(double *);
__declspec(dllimport) double __cdecl _Sinh(double, double);

__declspec(dllimport) short __cdecl _Exp(double *, double, short);
extern __declspec(dllimport)  _Dconst _Denorm, _Hugeval, _Inf,
	_Nan, _Snan;

		
__declspec(dllimport) float __cdecl _FCosh(float, float);
__declspec(dllimport) short __cdecl _FDtest(float *);
__declspec(dllimport) float __cdecl _FSinh(float, float);

__declspec(dllimport) short __cdecl _FExp(float *, float, short);
extern __declspec(dllimport)  _Dconst _FDenorm, _FInf, _FNan, _FSnan;

		
__declspec(dllimport) long double __cdecl _LCosh(long double, long double);
__declspec(dllimport) short __cdecl _LDtest(long double *);
__declspec(dllimport) long double __cdecl _LSinh(long double, long double);

__declspec(dllimport) short __cdecl _LExp(long double *, long double, short);
extern __declspec(dllimport)  _Dconst _LDenorm, _LInf, _LNan, _LSnan;

 
 
}
 
 

 
 #pragma warning(pop)
 #pragma pack(pop)











#pragma once













#pragma once




__pragma(pack(push, 8)) extern "C" {




    


        


            
        
    




















    






        
    



























































































































































    













__declspec(dllimport) unsigned int __cdecl _clearfp(void);

#pragma warning(push)
#pragma warning(disable: 4141)

 __declspec(deprecated("This function or variable may be unsafe. Consider using " "_controlfp_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__declspec(dllimport) unsigned int __cdecl _controlfp(
      unsigned int _NewValue,
      unsigned int _Mask
    );

#pragma warning(pop)


__declspec(dllimport) void __cdecl _set_controlfp(
      unsigned int _NewValue,
      unsigned int _Mask
    );


__declspec(dllimport) errno_t __cdecl _controlfp_s(
      unsigned int* _CurrentState,
           unsigned int  _NewValue,
           unsigned int  _Mask
    );


__declspec(dllimport) unsigned int __cdecl _statusfp(void);


__declspec(dllimport) void __cdecl _fpreset(void);















__declspec(dllimport) unsigned int __cdecl _control87(
      unsigned int _NewValue,
      unsigned int _Mask
    );












 
__declspec(dllimport) int* __cdecl __fpecode(void);



 
__declspec(dllimport) int __cdecl __fpe_flt_rounds(void);












  __declspec(dllimport) double __cdecl _copysign(  double _Number,   double _Sign);
  __declspec(dllimport) double __cdecl _chgsign(  double _X);
  __declspec(dllimport) double __cdecl _scalb(  double _X,   long _Y);
  __declspec(dllimport) double __cdecl _logb(  double _X);
  __declspec(dllimport) double __cdecl _nextafter(  double _X,   double _Y);
  __declspec(dllimport) int    __cdecl _finite(  double _X);
  __declspec(dllimport) int    __cdecl _isnan(  double _X);
  __declspec(dllimport) int    __cdecl _fpclass(  double _X);


      __declspec(dllimport) float __cdecl _scalbf(  float _X,   long _Y);











    
    
    

    
    __declspec(dllimport) void __cdecl fpreset(void);

    
    

    
    

    
    
    

    
    
    
    
    
    
    

    
    
    

    
    
    
    
    

    
    
    
    

    

    
    
    
    
    
    

    
    
    
    

    
    
    
    
    
    

    
    
    
    

    





} __pragma(pack(pop))










#pragma once














#pragma once









 












__pragma(pack(push, 8)) extern "C" {



#pragma warning(push)
#pragma warning(disable:4738) 
#pragma warning(disable:4820) 




    
    
    struct _exception
    {
        int    type;   
        char*  name;   
        double arg1;   
        double arg2;   
        double retval; 
    };

    
    
    
        

        struct _complex
        {
            double x, y; 
        };

        



    


typedef float float_t;
typedef double double_t;













    
        extern double const _HUGE;
    





    



















































































void __cdecl _fperrraise(  int _Except);

  __declspec(dllimport) short __cdecl _dclass(  double _X);
  __declspec(dllimport) short __cdecl _ldclass(  long double _X);
  __declspec(dllimport) short __cdecl _fdclass(  float _X);

  __declspec(dllimport) int __cdecl _dsign(  double _X);
  __declspec(dllimport) int __cdecl _ldsign(  long double _X);
  __declspec(dllimport) int __cdecl _fdsign(  float _X);

  __declspec(dllimport) int __cdecl _dpcomp(  double _X,   double _Y);
  __declspec(dllimport) int __cdecl _ldpcomp(  long double _X,   long double _Y);
  __declspec(dllimport) int __cdecl _fdpcomp(  float _X,   float _Y);

  __declspec(dllimport) short __cdecl _dtest(  double* _Px);
  __declspec(dllimport) short __cdecl _ldtest(  long double* _Px);
  __declspec(dllimport) short __cdecl _fdtest(  float* _Px);

__declspec(dllimport) short __cdecl _d_int(  double* _Px,   short _Xexp);
__declspec(dllimport) short __cdecl _ld_int(  long double* _Px,   short _Xexp);
__declspec(dllimport) short __cdecl _fd_int(  float* _Px,   short _Xexp);

__declspec(dllimport) short __cdecl _dscale(  double* _Px,   long _Lexp);
__declspec(dllimport) short __cdecl _ldscale(  long double* _Px,   long _Lexp);
__declspec(dllimport) short __cdecl _fdscale(  float* _Px,   long _Lexp);

__declspec(dllimport) short __cdecl _dunscale(  short* _Pex,   double* _Px);
__declspec(dllimport) short __cdecl _ldunscale(  short* _Pex,   long double* _Px);
__declspec(dllimport) short __cdecl _fdunscale(  short* _Pex,   float* _Px);

  __declspec(dllimport) short __cdecl _dexp(  double* _Px,   double _Y,   long _Eoff);
  __declspec(dllimport) short __cdecl _ldexp(  long double* _Px,   long double _Y,   long _Eoff);
  __declspec(dllimport) short __cdecl _fdexp(  float* _Px,   float _Y,   long _Eoff);

  __declspec(dllimport) short __cdecl _dnorm(  unsigned short* _Ps);
  __declspec(dllimport) short __cdecl _fdnorm(  unsigned short* _Ps);

  __declspec(dllimport) double __cdecl _dpoly(  double _X,   double const* _Tab,   int _N);
  __declspec(dllimport) long double __cdecl _ldpoly(  long double _X,   long double const* _Tab,   int _N);
  __declspec(dllimport) float __cdecl _fdpoly(  float _X,   float const* _Tab,   int _N);

  __declspec(dllimport) double __cdecl _dlog(  double _X,   int _Baseflag);
  __declspec(dllimport) long double __cdecl _ldlog(  long double _X,   int _Baseflag);
  __declspec(dllimport) float __cdecl _fdlog(  float _X,   int _Baseflag);

  __declspec(dllimport) double __cdecl _dsin(  double _X,   unsigned int _Qoff);
  __declspec(dllimport) long double __cdecl _ldsin(  long double _X,   unsigned int _Qoff);
  __declspec(dllimport) float __cdecl _fdsin(  float _X,   unsigned int _Qoff);


typedef union
{   
    unsigned short _Sh[4];
    double _Val;
} _double_val;


typedef union
{   
    unsigned short _Sh[2];
    float _Val;
} _float_val;


typedef union
{   
    unsigned short _Sh[4];
    long double _Val;
} _ldouble_val;

typedef union
{   
    unsigned short _Word[4];
    float _Float;
    double _Double;
    long double _Long_double;
} _float_const;

extern const _float_const _Denorm_C,  _Inf_C,  _Nan_C,  _Snan_C, _Hugeval_C;
extern const _float_const _FDenorm_C, _FInf_C, _FNan_C, _FSnan_C;
extern const _float_const _LDenorm_C, _LInf_C, _LNan_C, _LSnan_C;

extern const _float_const _Eps_C,  _Rteps_C;
extern const _float_const _FEps_C, _FRteps_C;
extern const _float_const _LEps_C, _LRteps_C;

extern const double      _Zero_C,  _Xbig_C;
extern const float       _FZero_C, _FXbig_C;
extern const long double _LZero_C, _LXbig_C;




























extern "C++"
{
      inline int fpclassify(  float _X) throw()
    {
        return _fdtest(&_X);
    }

      inline int fpclassify(  double _X) throw()
    {
        return _dtest(&_X);
    }

      inline int fpclassify(  long double _X) throw()
    {
        return _ldtest(&_X);
    }

      inline bool signbit(  float _X) throw()
    {
        return _fdsign(_X) != 0;
    }

      inline bool signbit(  double _X) throw()
    {
        return _dsign(_X) != 0;
    }

      inline bool signbit(  long double _X) throw()
    {
        return _ldsign(_X) != 0;
    }

      inline int _fpcomp(  float _X,   float _Y) throw()
    {
        return _fdpcomp(_X, _Y);
    }

      inline int _fpcomp(  double _X,   double _Y) throw()
    {
        return _dpcomp(_X, _Y);
    }

      inline int _fpcomp(  long double _X,   long double _Y) throw()
    {
        return _ldpcomp(_X, _Y);
    }

    template <class _Trc, class _Tre> struct _Combined_type
    {   
        typedef float _Type;
    };

    template <> struct _Combined_type<float, double>
    {   
        typedef double _Type;
    };

    template <> struct _Combined_type<float, long double>
    {   
        typedef long double _Type;
    };

    template <class _Ty, class _T2> struct _Real_widened
    {   
        typedef long double _Type;
    };

    template <> struct _Real_widened<float, float>
    {   
        typedef float _Type;
    };

    template <> struct _Real_widened<float, double>
    {   
        typedef double _Type;
    };

    template <> struct _Real_widened<double, float>
    {   
        typedef double _Type;
    };

    template <> struct _Real_widened<double, double>
    {   
        typedef double _Type;
    };

    template <class _Ty> struct _Real_type
    {   
        typedef double _Type;   
    };

    template <> struct _Real_type<float>
    {   
        typedef float _Type;
    };

    template <> struct _Real_type<long double>
    {   
        typedef long double _Type;
    };

    template <class _T1, class _T2>
      inline int _fpcomp(  _T1 _X,   _T2 _Y) throw()
    {   
        typedef typename _Combined_type<float,
            typename _Real_widened<
            typename _Real_type<_T1>::_Type,
            typename _Real_type<_T2>::_Type>::_Type>::_Type _Tw;
        return _fpcomp((_Tw)_X, (_Tw)_Y);
    }

    template <class _Ty>
      inline bool isfinite(  _Ty _X) throw()
    {
        return fpclassify(_X) <= 0;
    }

    template <class _Ty>
      inline bool isinf(  _Ty _X) throw()
    {
        return fpclassify(_X) == 1;
    }

    template <class _Ty>
      inline bool isnan(  _Ty _X) throw()
    {
        return fpclassify(_X) == 2;
    }

    template <class _Ty>
      inline bool isnormal(  _Ty _X) throw()
    {
        return fpclassify(_X) == (-1);
    }

    template <class _Ty1, class _Ty2>
      inline bool isgreater(  _Ty1 _X,   _Ty2 _Y) throw()
    {
        return (_fpcomp(_X, _Y) & 4) != 0;
    }

    template <class _Ty1, class _Ty2>
      inline bool isgreaterequal(  _Ty1 _X,   _Ty2 _Y) throw()
    {
        return (_fpcomp(_X, _Y) & (2 | 4)) != 0;
    }

    template <class _Ty1, class _Ty2>
      inline bool isless(  _Ty1 _X,   _Ty2 _Y) throw()
    {
        return (_fpcomp(_X, _Y) & 1) != 0;
    }

    template <class _Ty1, class _Ty2>
      inline bool islessequal(  _Ty1 _X,   _Ty2 _Y) throw()
    {
        return (_fpcomp(_X, _Y) & (1 | 2)) != 0;
    }

    template <class _Ty1, class _Ty2>
      inline bool islessgreater(  _Ty1 _X,   _Ty2 _Y) throw()
    {
        return (_fpcomp(_X, _Y) & (1 | 4)) != 0;
    }

    template <class _Ty1, class _Ty2>
      inline bool isunordered(  _Ty1 _X,   _Ty2 _Y) throw()
    {
        return _fpcomp(_X, _Y) == 0;
    }
}  






  int       __cdecl abs(  int _X);
  long      __cdecl labs(  long _X);
  long long __cdecl llabs(  long long _X);

  double __cdecl acos(  double _X);
  double __cdecl asin(  double _X);
  double __cdecl atan(  double _X);
  double __cdecl atan2(  double _Y,   double _X);

  double __cdecl cos(  double _X);
  double __cdecl cosh(  double _X);
  double __cdecl exp(  double _X);
   double __cdecl fabs(  double _X);
  double __cdecl fmod(  double _X,   double _Y);
  double __cdecl log(  double _X);
  double __cdecl log10(  double _X);
  double __cdecl pow(  double _X,   double _Y);
  double __cdecl sin(  double _X);
  double __cdecl sinh(  double _X);
   double __cdecl sqrt(  double _X);
  double __cdecl tan(  double _X);
  double __cdecl tanh(  double _X);

  __declspec(dllimport) double    __cdecl acosh(  double _X);
  __declspec(dllimport) double    __cdecl asinh(  double _X);
  __declspec(dllimport) double    __cdecl atanh(  double _X);
  __declspec(dllimport)  double    __cdecl atof(  char const* _String);
  __declspec(dllimport)  double    __cdecl _atof_l(  char const* _String,   _locale_t _Locale);
  __declspec(dllimport) double    __cdecl _cabs(  struct _complex _Complex_value);
  __declspec(dllimport) double    __cdecl cbrt(  double _X);
  __declspec(dllimport) double    __cdecl ceil(  double _X);
  __declspec(dllimport) double    __cdecl _chgsign(  double _X);
  __declspec(dllimport) double    __cdecl copysign(  double _Number,   double _Sign);
  __declspec(dllimport) double    __cdecl _copysign(  double _Number,   double _Sign);
  __declspec(dllimport) double    __cdecl erf(  double _X);
  __declspec(dllimport) double    __cdecl erfc(  double _X);
  __declspec(dllimport) double    __cdecl exp2(  double _X);
  __declspec(dllimport) double    __cdecl expm1(  double _X);
  __declspec(dllimport) double    __cdecl fdim(  double _X,   double _Y);
  __declspec(dllimport) double    __cdecl floor(  double _X);
  __declspec(dllimport) double    __cdecl fma(  double _X,   double _Y,   double _Z);
  __declspec(dllimport) double    __cdecl fmax(  double _X,   double _Y);
  __declspec(dllimport) double    __cdecl fmin(  double _X,   double _Y);
  __declspec(dllimport) double    __cdecl frexp(  double _X,   int* _Y);
  __declspec(dllimport) double    __cdecl hypot(  double _X,   double _Y);
  __declspec(dllimport) double    __cdecl _hypot(  double _X,   double _Y);
  __declspec(dllimport) int       __cdecl ilogb(  double _X);
  __declspec(dllimport) double    __cdecl ldexp(  double _X,   int _Y);
  __declspec(dllimport) double    __cdecl lgamma(  double _X);
  __declspec(dllimport) long long __cdecl llrint(  double _X);
  __declspec(dllimport) long long __cdecl llround(  double _X);
  __declspec(dllimport) double    __cdecl log1p(  double _X);
  __declspec(dllimport) double    __cdecl log2(  double _X);
  __declspec(dllimport) double    __cdecl logb(  double _X);
  __declspec(dllimport) long      __cdecl lrint(  double _X);
  __declspec(dllimport) long      __cdecl lround(  double _X);

int __cdecl _matherr(  struct _exception* _Except);

  __declspec(dllimport) double __cdecl modf(  double _X,   double* _Y);
  __declspec(dllimport) double __cdecl nan(  char const*);
  __declspec(dllimport) double __cdecl nearbyint(  double _X);
  __declspec(dllimport) double __cdecl nextafter(  double _X,   double _Y);
  __declspec(dllimport) double __cdecl nexttoward(  double _X,   long double _Y);
  __declspec(dllimport) double __cdecl remainder(  double _X,   double _Y);
  __declspec(dllimport) double __cdecl remquo(  double _X,   double _Y,   int* _Z);
  __declspec(dllimport) double __cdecl rint(  double _X);
  __declspec(dllimport) double __cdecl round(  double _X);
  __declspec(dllimport) double __cdecl scalbln(  double _X,   long _Y);
  __declspec(dllimport) double __cdecl scalbn(  double _X,   int _Y);
  __declspec(dllimport) double __cdecl tgamma(  double _X);
  __declspec(dllimport) double __cdecl trunc(  double _X);
  __declspec(dllimport) double __cdecl _j0(  double _X );
  __declspec(dllimport) double __cdecl _j1(  double _X );
  __declspec(dllimport) double __cdecl _jn(int _X,   double _Y);
  __declspec(dllimport) double __cdecl _y0(  double _X);
  __declspec(dllimport) double __cdecl _y1(  double _X);
  __declspec(dllimport) double __cdecl _yn(  int _X,   double _Y);

  __declspec(dllimport) float     __cdecl acoshf(  float _X);
  __declspec(dllimport) float     __cdecl asinhf(  float _X);
  __declspec(dllimport) float     __cdecl atanhf(  float _X);
  __declspec(dllimport) float     __cdecl cbrtf(  float _X);
  __declspec(dllimport) float     __cdecl _chgsignf(  float _X);
  __declspec(dllimport) float     __cdecl copysignf(  float _Number,   float _Sign);
  __declspec(dllimport) float     __cdecl _copysignf(  float _Number,   float _Sign);
  __declspec(dllimport) float     __cdecl erff(  float _X);
  __declspec(dllimport) float     __cdecl erfcf(  float _X);
  __declspec(dllimport) float     __cdecl expm1f(  float _X);
  __declspec(dllimport) float     __cdecl exp2f(  float _X);
  __declspec(dllimport) float     __cdecl fdimf(  float _X,   float _Y);
  __declspec(dllimport) float     __cdecl fmaf(  float _X,   float _Y,   float _Z);
  __declspec(dllimport) float     __cdecl fmaxf(  float _X,   float _Y);
  __declspec(dllimport) float     __cdecl fminf(  float _X,   float _Y);
  __declspec(dllimport) float     __cdecl _hypotf(  float _X,   float _Y);
  __declspec(dllimport) int       __cdecl ilogbf(  float _X);
  __declspec(dllimport) float     __cdecl lgammaf(  float _X);
  __declspec(dllimport) long long __cdecl llrintf(  float _X);
  __declspec(dllimport) long long __cdecl llroundf(  float _X);
  __declspec(dllimport) float     __cdecl log1pf(  float _X);
  __declspec(dllimport) float     __cdecl log2f(  float _X);
  __declspec(dllimport) float     __cdecl logbf(  float _X);
  __declspec(dllimport) long      __cdecl lrintf(  float _X);
  __declspec(dllimport) long      __cdecl lroundf(  float _X);
  __declspec(dllimport) float     __cdecl nanf(  char const*);
  __declspec(dllimport) float     __cdecl nearbyintf(  float _X);
  __declspec(dllimport) float     __cdecl nextafterf(  float _X,   float _Y);
  __declspec(dllimport) float     __cdecl nexttowardf(  float _X,   long double _Y);
  __declspec(dllimport) float     __cdecl remainderf(  float _X,   float _Y);
  __declspec(dllimport) float     __cdecl remquof(  float _X,   float _Y,   int* _Z);
  __declspec(dllimport) float     __cdecl rintf(  float _X);
  __declspec(dllimport) float     __cdecl roundf(  float _X);
  __declspec(dllimport) float     __cdecl scalblnf(  float _X,   long _Y);
  __declspec(dllimport) float     __cdecl scalbnf(  float _X,   int _Y);
  __declspec(dllimport) float     __cdecl tgammaf(  float _X);
  __declspec(dllimport) float     __cdecl truncf(  float _X);







      __declspec(dllimport) float __cdecl _logbf(  float _X);
      __declspec(dllimport) float __cdecl _nextafterf(  float _X,   float _Y);
      __declspec(dllimport) int   __cdecl _finitef(  float _X);
      __declspec(dllimport) int   __cdecl _isnanf(  float _X);
      __declspec(dllimport) int   __cdecl _fpclassf(  float _X);

      __declspec(dllimport) int   __cdecl _set_FMA3_enable(  int _Flag);
      __declspec(dllimport) int   __cdecl _get_FMA3_enable(void);












      __declspec(dllimport) float __cdecl acosf(  float _X);
      __declspec(dllimport) float __cdecl asinf(  float _X);
      __declspec(dllimport) float __cdecl atan2f(  float _Y,   float _X);
      __declspec(dllimport) float __cdecl atanf(  float _X);
      __declspec(dllimport) float __cdecl ceilf(  float _X);
      __declspec(dllimport) float __cdecl cosf(  float _X);
      __declspec(dllimport) float __cdecl coshf(  float _X);
      __declspec(dllimport) float __cdecl expf(  float _X);



















































      __inline float __cdecl fabsf(  float _X)
    {
        return (float)fabs(_X);
    }





      __declspec(dllimport) float __cdecl floorf(  float _X);
      __declspec(dllimport) float __cdecl fmodf(  float _X,   float _Y);















  __inline float __cdecl frexpf(  float _X,   int *_Y)
{
    return (float)frexp(_X, _Y);
}

  __inline float __cdecl hypotf(  float _X,   float _Y)
{
    return _hypotf(_X, _Y);
}

  __inline float __cdecl ldexpf(  float _X,   int _Y)
{
    return (float)ldexp(_X, _Y);
}



      __declspec(dllimport) float  __cdecl log10f(  float _X);
      __declspec(dllimport) float  __cdecl logf(  float _X);
      __declspec(dllimport) float  __cdecl modff(  float _X,   float *_Y);
      __declspec(dllimport) float  __cdecl powf(  float _X,   float _Y);
      __declspec(dllimport) float  __cdecl sinf(  float _X);
      __declspec(dllimport) float  __cdecl sinhf(  float _X);
      __declspec(dllimport) float  __cdecl sqrtf(  float _X);
      __declspec(dllimport) float  __cdecl tanf(  float _X);
      __declspec(dllimport) float  __cdecl tanhf(  float _X);





















































  __declspec(dllimport) long double __cdecl acoshl(  long double _X);

  __inline long double __cdecl acosl(  long double _X)
{
    return acos((double)_X);
}

  __declspec(dllimport) long double __cdecl asinhl(  long double _X);

  __inline long double __cdecl asinl(  long double _X)
{
    return asin((double)_X);
}

  __inline long double __cdecl atan2l(  long double _Y,   long double _X)
{
    return atan2((double)_Y, (double)_X);
}

  __declspec(dllimport) long double __cdecl atanhl(  long double _X);

  __inline long double __cdecl atanl(  long double _X)
{
    return atan((double)_X);
}

  __declspec(dllimport) long double __cdecl cbrtl(  long double _X);

  __inline long double __cdecl ceill(  long double _X)
{
    return ceil((double)_X);
}

  __inline long double __cdecl _chgsignl(  long double _X)
{
    return _chgsign((double)_X);
}

  __declspec(dllimport) long double __cdecl copysignl(  long double _Number,   long double _Sign);

  __inline long double __cdecl _copysignl(  long double _Number,   long double _Sign)
{
    return _copysign((double)_Number, (double)_Sign);
}

  __inline long double __cdecl coshl(  long double _X)
{
    return cosh((double)_X);
}

  __inline long double __cdecl cosl(  long double _X)
{
    return cos((double)_X);
}

  __declspec(dllimport) long double __cdecl erfl(  long double _X);
  __declspec(dllimport) long double __cdecl erfcl(  long double _X);

  __inline long double __cdecl expl(  long double _X)
{
    return exp((double)_X);
}

  __declspec(dllimport) long double __cdecl exp2l(  long double _X);
  __declspec(dllimport) long double __cdecl expm1l(  long double _X);

  __inline long double __cdecl fabsl(  long double _X)
{
    return fabs((double)_X);
}

  __declspec(dllimport) long double __cdecl fdiml(  long double _X,   long double _Y);

  __inline long double __cdecl floorl(  long double _X)
{
    return floor((double)_X);
}

  __declspec(dllimport) long double __cdecl fmal(  long double _X,   long double _Y,   long double _Z);
  __declspec(dllimport) long double __cdecl fmaxl(  long double _X,   long double _Y);
  __declspec(dllimport) long double __cdecl fminl(  long double _X,   long double _Y);

  __inline long double __cdecl fmodl(  long double _X,   long double _Y)
{
    return fmod((double)_X, (double)_Y);
}

  __inline long double __cdecl frexpl(  long double _X,   int *_Y)
{
    return frexp((double)_X, _Y);
}

  __declspec(dllimport) int __cdecl ilogbl(  long double _X);

  __inline long double __cdecl _hypotl(  long double _X,   long double _Y)
{
    return _hypot((double)_X, (double)_Y);
}

  __inline long double __cdecl hypotl(  long double _X,   long double _Y)
{
    return _hypot((double)_X, (double)_Y);
}

  __inline long double __cdecl ldexpl(  long double _X,   int _Y)
{
    return ldexp((double)_X, _Y);
}

  __declspec(dllimport) long double __cdecl lgammal(  long double _X);
  __declspec(dllimport) long long __cdecl llrintl(  long double _X);
  __declspec(dllimport) long long __cdecl llroundl(  long double _X);

  __inline long double __cdecl logl(  long double _X)
{
    return log((double)_X);
}

  __inline long double __cdecl log10l(  long double _X)
{
    return log10((double)_X);
}

  __declspec(dllimport) long double __cdecl log1pl(  long double _X);
  __declspec(dllimport) long double __cdecl log2l(  long double _X);
  __declspec(dllimport) long double __cdecl logbl(  long double _X);
  __declspec(dllimport) long __cdecl lrintl(  long double _X);
  __declspec(dllimport) long __cdecl lroundl(  long double _X);

  __inline long double __cdecl modfl(  long double _X,   long double* _Y)
{
    double _F, _I;
    _F = modf((double)_X, &_I);
    *_Y = _I;
    return _F;
}

  __declspec(dllimport) long double __cdecl nanl(  char const*);
  __declspec(dllimport) long double __cdecl nearbyintl(  long double _X);
  __declspec(dllimport) long double __cdecl nextafterl(  long double _X,   long double _Y);
  __declspec(dllimport) long double __cdecl nexttowardl(  long double _X,   long double _Y);

  __inline long double __cdecl powl(  long double _X,   long double _Y)
{
    return pow((double)_X, (double)_Y);
}

  __declspec(dllimport) long double __cdecl remainderl(  long double _X,   long double _Y);
  __declspec(dllimport) long double __cdecl remquol(  long double _X,   long double _Y,   int* _Z);
  __declspec(dllimport) long double __cdecl rintl(  long double _X);
  __declspec(dllimport) long double __cdecl roundl(  long double _X);
  __declspec(dllimport) long double __cdecl scalblnl(  long double _X,   long _Y);
  __declspec(dllimport) long double __cdecl scalbnl(  long double _X,   int _Y);

  __inline long double __cdecl sinhl(  long double _X)
{
    return sinh((double)_X);
}

  __inline long double __cdecl sinl(  long double _X)
{
    return sin((double)_X);
}

  __inline long double __cdecl sqrtl(  long double _X)
{
    return sqrt((double)_X);
}

  __inline long double __cdecl tanhl(  long double _X)
{
    return tanh((double)_X);
}

  __inline long double __cdecl tanl(  long double _X)
{
    return tan((double)_X);
}

  __declspec(dllimport) long double __cdecl tgammal(  long double _X);
  __declspec(dllimport) long double __cdecl truncl(  long double _X);









    
    
    
    
    
    
    

    

    
        
            extern double HUGE;
        



        __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_j0" ". See online help for details."))   __declspec(dllimport) double __cdecl j0(  double _X);
        __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_j1" ". See online help for details."))   __declspec(dllimport) double __cdecl j1(  double _X);
        __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_jn" ". See online help for details."))   __declspec(dllimport) double __cdecl jn(  int _X,   double _Y);
        __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_y0" ". See online help for details."))   __declspec(dllimport) double __cdecl y0(  double _X);
        __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_y1" ". See online help for details."))   __declspec(dllimport) double __cdecl y1(  double _X);
        __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_yn" ". See online help for details."))   __declspec(dllimport) double __cdecl yn(  int _X,   double _Y);
    




#pragma warning(pop)



} __pragma(pack(pop))






    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    





 

  inline double abs(  double _Xx) noexcept
	{
	return (:: fabs(_Xx));
	}

  inline double pow(  double _Xx,   int _Yx) noexcept
	{
	if (_Yx == 2)
		return (_Xx * _Xx);

	return (:: pow(_Xx, static_cast<double>(_Yx)));
	}

  inline float abs(  float _Xx) noexcept
	{
	return (:: fabsf(_Xx));
	}

  inline float acos(  float _Xx) noexcept
	{
	return (:: acosf(_Xx));
	}

  inline float acosh(  float _Xx) noexcept
	{
	return (:: acoshf(_Xx));
	}

  inline float asin(  float _Xx) noexcept
	{
	return (:: asinf(_Xx));
	}

  inline float asinh(  float _Xx) noexcept
	{
	return (:: asinhf(_Xx));
	}

  inline float atan(  float _Xx) noexcept
	{
	return (:: atanf(_Xx));
	}

  inline float atanh(  float _Xx) noexcept
	{
	return (:: atanhf(_Xx));
	}

  inline float atan2(  float _Yx,   float _Xx) noexcept
	{
	return (:: atan2f(_Yx, _Xx));
	}

  inline float cbrt(  float _Xx) noexcept
	{
	return (:: cbrtf(_Xx));
	}

  inline float ceil(  float _Xx) noexcept
	{
	return (:: ceilf(_Xx));
	}

  inline float copysign(  float _Number,
	  float _Sign) noexcept
	{
	return (:: copysignf(_Number, _Sign));
	}

  inline float cos(  float _Xx) noexcept
	{
	return (:: cosf(_Xx));
	}

  inline float cosh(  float _Xx) noexcept
	{
	return (:: coshf(_Xx));
	}

  inline float erf(  float _Xx) noexcept
	{
	return (:: erff(_Xx));
	}

  inline float erfc(  float _Xx) noexcept
	{
	return (:: erfcf(_Xx));
	}

  inline float exp(  float _Xx) noexcept
	{
	return (:: expf(_Xx));
	}

  inline float exp2(  float _Xx) noexcept
	{
	return (:: exp2f(_Xx));
	}

  inline float expm1(  float _Xx) noexcept
	{
	return (:: expm1f(_Xx));
	}

  inline float fabs(  float _Xx) noexcept
	{
	return (:: fabsf(_Xx));
	}

  inline float fdim(  float _Xx,   float _Yx) noexcept
	{
	return (:: fdimf(_Xx, _Yx));
	}

  inline float floor(  float _Xx) noexcept
	{
	return (:: floorf(_Xx));
	}

  inline float fma(  float _Xx,   float _Yx,
	  float _Zx) noexcept
	{
	return (:: fmaf(_Xx, _Yx, _Zx));
	}

  inline float fmax(  float _Xx,   float _Yx) noexcept
	{
	return (:: fmaxf(_Xx, _Yx));
	}

  inline float fmin(  float _Xx,   float _Yx) noexcept
	{
	return (:: fminf(_Xx, _Yx));
	}

  inline float fmod(  float _Xx,   float _Yx) noexcept
	{
	return (:: fmodf(_Xx, _Yx));
	}

  inline float frexp(  float _Xx,   int* _Yx) noexcept
	{
	return (:: frexpf(_Xx, _Yx));
	}

  inline float hypot(  float _Xx,   float _Yx) noexcept
	{
	return (:: hypotf(_Xx, _Yx));
	}

  inline int ilogb(  float _Xx) noexcept
	{
	return (:: ilogbf(_Xx));
	}

  inline float ldexp(  float _Xx,   int _Yx) noexcept
	{
	return (:: ldexpf(_Xx, _Yx));
	}

  inline float lgamma(  float _Xx) noexcept
	{
	return (:: lgammaf(_Xx));
	}

  inline long long llrint(  float _Xx) noexcept
	{
	return (:: llrintf(_Xx));
	}

  inline long long llround(  float _Xx) noexcept
	{
	return (:: llroundf(_Xx));
	}

  inline float log(  float _Xx) noexcept
	{
	return (:: logf(_Xx));
	}

  inline float log10(  float _Xx) noexcept
	{
	return (:: log10f(_Xx));
	}

  inline float log1p(  float _Xx) noexcept
	{
	return (:: log1pf(_Xx));
	}

  inline float log2(  float _Xx) noexcept
	{
	return (:: log2f(_Xx));
	}

  inline float logb(  float _Xx) noexcept
	{
	return (:: logbf(_Xx));
	}

  inline long lrint(  float _Xx) noexcept
	{
	return (:: lrintf(_Xx));
	}

  inline long lround(  float _Xx) noexcept
	{
	return (:: lroundf(_Xx));
	}

  inline float modf(  float _Xx,   float* _Yx) noexcept
	{
	return (:: modff(_Xx, _Yx));
	}

  inline float nearbyint(  float _Xx) noexcept
	{
	return (:: nearbyintf(_Xx));
	}

  inline float nextafter(  float _Xx,   float _Yx) noexcept
	{
	return (:: nextafterf(_Xx, _Yx));
	}

  inline float nexttoward(  float _Xx,
	  long double _Yx) noexcept
	{
	return (:: nexttowardf(_Xx, _Yx));
	}

  inline float pow(  float _Xx,
	  float _Yx) noexcept
	{
	return (:: powf(_Xx, _Yx));
	}

  inline float pow(  float _Xx,   int _Yx) noexcept
	{
	if (_Yx == 2)
		return (_Xx * _Xx);

	return (:: powf(_Xx, static_cast<float>(_Yx)));
	}

  inline float remainder(  float _Xx,   float _Yx) noexcept
	{
	return (:: remainderf(_Xx, _Yx));
	}

  inline float remquo(  float _Xx,   float _Yx,
	  int *_Zx) noexcept
	{
	return (:: remquof(_Xx, _Yx, _Zx));
	}

  inline float rint(  float _Xx) noexcept
	{
	return (:: rintf(_Xx));
	}

  inline float round(  float _Xx) noexcept
	{
	return (:: roundf(_Xx));
	}

  inline float scalbln(  float _Xx,   long _Yx) noexcept
	{
	return (:: scalblnf(_Xx, _Yx));
	}

  inline float scalbn(  float _Xx,   int _Yx) noexcept
	{
	return (:: scalbnf(_Xx, _Yx));
	}

  inline float sin(  float _Xx) noexcept
	{
	return (:: sinf(_Xx));
	}

  inline float sinh(  float _Xx) noexcept
	{
	return (:: sinhf(_Xx));
	}

  inline float sqrt(  float _Xx) noexcept
	{
	return (:: sqrtf(_Xx));
	}

  inline float tan(  float _Xx) noexcept
	{
	return (:: tanf(_Xx));
	}

  inline float tanh(  float _Xx) noexcept
	{
	return (:: tanhf(_Xx));
	}

  inline float tgamma(  float _Xx) noexcept
	{
	return (:: tgammaf(_Xx));
	}

  inline float trunc(  float _Xx) noexcept
	{
	return (:: truncf(_Xx));
	}

  inline long double abs(  long double _Xx) noexcept
	{
	return (:: fabsl(_Xx));
	}

  inline long double acos(  long double _Xx) noexcept
	{
	return (:: acosl(_Xx));
	}

  inline long double acosh(  long double _Xx) noexcept
	{
	return (:: acoshl(_Xx));
	}

  inline long double asin(  long double _Xx) noexcept
	{
	return (:: asinl(_Xx));
	}

  inline long double asinh(  long double _Xx) noexcept
	{
	return (:: asinhl(_Xx));
	}

  inline long double atan(  long double _Xx) noexcept
	{
	return (:: atanl(_Xx));
	}

  inline long double atanh(  long double _Xx) noexcept
	{
	return (:: atanhl(_Xx));
	}

  inline long double atan2(  long double _Yx,
	  long double _Xx) noexcept
	{
	return (:: atan2l(_Yx, _Xx));
	}

  inline long double cbrt(  long double _Xx) noexcept
	{
	return (:: cbrtl(_Xx));
	}

  inline long double ceil(  long double _Xx) noexcept
	{
	return (:: ceill(_Xx));
	}

  inline long double copysign(  long double _Number,
	  long double _Sign) noexcept
	{
	return (:: copysignl(_Number, _Sign));
	}

  inline long double cos(  long double _Xx) noexcept
	{
	return (:: cosl(_Xx));
	}

  inline long double cosh(  long double _Xx) noexcept
	{
	return (:: coshl(_Xx));
	}

  inline long double erf(  long double _Xx) noexcept
	{
	return (:: erfl(_Xx));
	}

  inline long double erfc(  long double _Xx) noexcept
	{
	return (:: erfcl(_Xx));
	}

  inline long double exp(  long double _Xx) noexcept
	{
	return (:: expl(_Xx));
	}

  inline long double exp2(  long double _Xx) noexcept
	{
	return (:: exp2l(_Xx));
	}

  inline long double expm1(  long double _Xx) noexcept
	{
	return (:: expm1l(_Xx));
	}

  inline long double fabs(  long double _Xx) noexcept
	{
	return (:: fabsl(_Xx));
	}

  inline long double fdim(  long double _Xx,
	  long double _Yx) noexcept
	{
	return (:: fdiml(_Xx, _Yx));
	}

  inline long double floor(  long double _Xx) noexcept
	{
	return (:: floorl(_Xx));
	}

  inline long double fma(  long double _Xx,
	  long double _Yx,   long double _Zx) noexcept
	{
	return (:: fmal(_Xx, _Yx, _Zx));
	}

  inline long double fmax(  long double _Xx,
	  long double _Yx) noexcept
	{
	return (:: fmaxl(_Xx, _Yx));
	}

  inline long double fmin(  long double _Xx,
	  long double _Yx) noexcept
	{
	return (:: fminl(_Xx, _Yx));
	}

  inline long double fmod(  long double _Xx,
	  long double _Yx) noexcept
	{
	return (:: fmodl(_Xx, _Yx));
	}

  inline long double frexp(  long double _Xx,
	  int* _Yx) noexcept
	{
	return (:: frexpl(_Xx, _Yx));
	}

  inline long double hypot(  long double _Xx,
	  long double _Yx) noexcept
	{
	return (:: hypotl(_Xx, _Yx));
	}

  inline int ilogb(  long double _Xx) noexcept
	{
	return (:: ilogbl(_Xx));
	}

  inline long double ldexp(  long double _Xx,
	  int _Yx) noexcept
	{
	return (:: ldexpl(_Xx, _Yx));
	}

  inline long double lgamma(  long double _Xx) noexcept
	{
	return (:: lgammal(_Xx));
	}

  inline long long llrint(  long double _Xx) noexcept
	{
	return (:: llrintl(_Xx));
	}

  inline long long llround(  long double _Xx) noexcept
	{
	return (:: llroundl(_Xx));
	}

  inline long double log(  long double _Xx) noexcept
	{
	return (:: logl(_Xx));
	}

  inline long double log10(  long double _Xx) noexcept
	{
	return (:: log10l(_Xx));
	}

  inline long double log1p(  long double _Xx) noexcept
	{
	return (:: log1pl(_Xx));
	}

  inline long double log2(  long double _Xx) noexcept
	{
	return (:: log2l(_Xx));
	}

  inline long double logb(  long double _Xx) noexcept
	{
	return (:: logbl(_Xx));
	}

  inline long lrint(  long double _Xx) noexcept
	{
	return (:: lrintl(_Xx));
	}

  inline long lround(  long double _Xx) noexcept
	{
	return (:: lroundl(_Xx));
	}

  inline long double modf(  long double _Xx,
	  long double* _Yx) noexcept
	{
	return (:: modfl(_Xx, _Yx));
	}

  inline long double nearbyint(  long double _Xx) noexcept
	{
	return (:: nearbyintl(_Xx));
	}

  inline long double nextafter(  long double _Xx,
	  long double _Yx) noexcept
	{
	return (:: nextafterl(_Xx, _Yx));
	}

  inline long double nexttoward(  long double _Xx,
	  long double _Yx) noexcept
	{
	return (:: nexttowardl(_Xx, _Yx));
	}

  inline long double pow(  long double _Xx,
	  long double _Yx) noexcept
	{
	return (:: powl(_Xx, _Yx));
	}

  inline long double pow(  long double _Xx,
	  int _Yx) noexcept
	{
	if (_Yx == 2)
		return (_Xx * _Xx);

	return (:: powl(_Xx, static_cast<long double>(_Yx)));
	}

  inline long double remainder(  long double _Xx,
	  long double _Yx) noexcept
	{
	return (:: remainderl(_Xx, _Yx));
	}

  inline long double remquo(  long double _Xx,
	  long double _Yx,   int *_Zx) noexcept
	{
	return (:: remquol(_Xx, _Yx, _Zx));
	}

  inline long double rint(  long double _Xx) noexcept
	{
	return (:: rintl(_Xx));
	}

  inline long double round(  long double _Xx) noexcept
	{
	return (:: roundl(_Xx));
	}

  inline long double scalbln(  long double _Xx,
	  long _Yx) noexcept
	{
	return (:: scalblnl(_Xx, _Yx));
	}

  inline long double scalbn(  long double _Xx,
	  int _Yx) noexcept
	{
	return (:: scalbnl(_Xx, _Yx));
	}

  inline long double sin(  long double _Xx) noexcept
	{
	return (:: sinl(_Xx));
	}

  inline long double sinh(  long double _Xx) noexcept
	{
	return (:: sinhl(_Xx));
	}

  inline long double sqrt(  long double _Xx) noexcept
	{
	return (:: sqrtl(_Xx));
	}

  inline long double tan(  long double _Xx) noexcept
	{
	return (:: tanl(_Xx));
	}

  inline long double tanh(  long double _Xx) noexcept
	{
	return (:: tanhl(_Xx));
	}

  inline long double tgamma(  long double _Xx) noexcept
	{
	return (:: tgammal(_Xx));
	}

  inline long double trunc(  long double _Xx) noexcept
	{
	return (:: truncl(_Xx));
	}

 






 


 
#pragma once





#pragma once





 #pragma pack(push,8)
 #pragma warning(push,3)
 
 

namespace std {
	
struct _Nil
	{	
	};

	
template<class _Ty,
	_Ty _Val>
	struct integral_constant
	{	
	static constexpr _Ty value = _Val;

	typedef _Ty value_type;
	typedef integral_constant<_Ty, _Val> type;

	constexpr operator value_type() const noexcept
		{	
		return (value);
		}

	constexpr value_type operator()() const noexcept
		{	
		return (value);
		}
	};

typedef integral_constant<bool, true> true_type;
typedef integral_constant<bool, false> false_type;

	
template<bool _Val>
	using bool_constant = integral_constant<bool, _Val>;

	
template<bool _Val>
	struct _Cat_base
		: integral_constant<bool, _Val>
	{	
	};

	
template<bool _Test,
	class _Ty = void>
	struct enable_if
	{	
	};

template<class _Ty>
	struct enable_if<true, _Ty>
	{	
	typedef _Ty type;
	};

	
template<bool _Test,
	class _Ty1,
	class _Ty2>
	struct conditional
	{	
	typedef _Ty2 type;
	};

template<class _Ty1,
	class _Ty2>
	struct conditional<true, _Ty1, _Ty2>
	{	
	typedef _Ty1 type;
	};

	
template<class _Ty1,
	class _Ty2>
	struct is_same
		: false_type
	{	
	};

template<class _Ty1>
	struct is_same<_Ty1, _Ty1>
		: true_type
	{	
	};

 
template<class _Ty,
	class _Uty>
	constexpr bool is_same_v = is_same<_Ty, _Uty>::value;
 

	
template<class _Ty>
	struct remove_const
	{	
	typedef _Ty type;
	};

template<class _Ty>
	struct remove_const<const _Ty>
	{	
	typedef _Ty type;
	};

	
template<class _Ty>
	struct remove_volatile
	{	
	typedef _Ty type;
	};

template<class _Ty>
	struct remove_volatile<volatile _Ty>
	{	
	typedef _Ty type;
	};

	
template<class _Ty>
	struct remove_cv
	{	
	typedef typename remove_const<typename remove_volatile<_Ty>::type>::type
		type;
	};

	
template<class _Ty>
	struct _Is_integral
		: false_type
	{	
	};

template<>
	struct _Is_integral<bool>
		: true_type
	{	
	};

template<>
	struct _Is_integral<char>
		: true_type
	{	
	};

template<>
	struct _Is_integral<unsigned char>
		: true_type
	{	
	};

template<>
	struct _Is_integral<signed char>
		: true_type
	{	
	};

 
template<>
	struct _Is_integral<wchar_t>
		: true_type
	{	
	};
 

template<>
	struct _Is_integral<unsigned short>
		: true_type
	{	
	};

template<>
	struct _Is_integral<signed short>
		: true_type
	{	
	};

template<>
	struct _Is_integral<unsigned int>
		: true_type
	{	
	};

template<>
	struct _Is_integral<signed int>
		: true_type
	{	
	};

template<>
	struct _Is_integral<unsigned long>
		: true_type
	{	
	};

template<>
	struct _Is_integral<signed long>
		: true_type
	{	
	};

template<>
	struct _Is_integral<char16_t>
		: true_type
	{	
	};

template<>
	struct _Is_integral<char32_t>
		: true_type
	{	
	};

template<>
	struct _Is_integral<long long>
		: true_type
	{	
	};

template<>
	struct _Is_integral<unsigned long long>
		: true_type
	{	
	};

	
template<class _Ty>
	struct is_integral
		: _Is_integral<typename remove_cv<_Ty>::type>
	{	
	};

 
template<class _Ty>
	constexpr bool is_integral_v = is_integral<_Ty>::value;
 

	
template<class _Ty>
	struct _Is_floating_point
		: false_type
	{	
	};

template<>
	struct _Is_floating_point<float>
		: true_type
	{	
	};

template<>
	struct _Is_floating_point<double>
		: true_type
	{	
	};

template<>
	struct _Is_floating_point<long double>
		: true_type
	{	
	};

	
template<class _Ty>
	struct is_floating_point
		: _Is_floating_point<typename remove_cv<_Ty>::type>
	{	
	};

 
template<class _Ty>
	constexpr bool is_floating_point_v = is_floating_point<_Ty>::value;
 

	
template<class _Ty>
	struct is_arithmetic
		: _Cat_base<is_integral<_Ty>::value
			|| is_floating_point<_Ty>::value>
	{	
	};

 
template<class _Ty>
	constexpr bool is_arithmetic_v = is_arithmetic<_Ty>::value;
 

	
template<class _Ty>
	struct remove_reference
	{	
	typedef _Ty type;
	};

template<class _Ty>
	struct remove_reference<_Ty&>
	{	
	typedef _Ty type;
	};

template<class _Ty>
	struct remove_reference<_Ty&&>
	{	
	typedef _Ty type;
	};

	
struct _Wrap_int
	{	
	_Wrap_int(int)
		{	
		}
	};

template<class _Ty>
	struct _Identity
	{	
	typedef _Ty type;
	};














}
 
 #pragma warning(pop)
 #pragma pack(pop)









 #pragma pack(push,8)
 #pragma warning(push,3)
 
 

namespace std {
template<class _Ty>
	struct _Promote_to_float
	{	
	typedef typename conditional<is_integral<_Ty>::value,
		double, _Ty>::type type;
	};

template<class _Ty1,
	class _Ty2>
	struct _Common_float_type
	{	
	typedef typename _Promote_to_float<_Ty1>::type _Ty1f;
	typedef typename _Promote_to_float<_Ty2>::type _Ty2f;
	typedef typename conditional<is_same<_Ty1f, long double>::value
		|| is_same<_Ty2f, long double>::value, long double,
		typename conditional<is_same<_Ty1f, double>::value
			|| is_same<_Ty2f, double>::value, double,
			float>::type>::type type;
	};
}








































template<class _Ty1,
	class _Ty2> inline
	typename ::std:: enable_if< ::std:: is_arithmetic<_Ty1>::value
		&& ::std:: is_arithmetic<_Ty2>::value,
		typename ::std:: _Common_float_type<_Ty1, _Ty2>::type>::type
	pow(const _Ty1 _Left, const _Ty2 _Right)
	{	
	typedef typename ::std:: _Common_float_type<_Ty1, _Ty2>::type type;
	return (:: pow(type(_Left), type(_Right)));
	}


extern "C"    double __cdecl acos(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type acos(_Ty _Left) { return (:: acos((double)_Left)); }
extern "C"    double __cdecl asin(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type asin(_Ty _Left) { return (:: asin((double)_Left)); }
extern "C"    double __cdecl atan(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type atan(_Ty _Left) { return (:: atan((double)_Left)); }
extern "C"    double __cdecl atan2(   double,   double); template<class _Ty1, class _Ty2> inline typename ::std:: enable_if< ::std:: is_arithmetic<_Ty1>::value && ::std:: is_arithmetic<_Ty2>::value, typename ::std:: _Common_float_type<_Ty1, _Ty2>::type>::type atan2(_Ty1 _Left, _Ty2 _Right) { typedef typename ::std:: _Common_float_type<_Ty1, _Ty2>::type type; return (:: atan2((type)_Left, (type)_Right)); }
extern "C"   __declspec(dllimport) double __cdecl ceil(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type ceil(_Ty _Left) { return (:: ceil((double)_Left)); }
extern "C"    double __cdecl cos(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type cos(_Ty _Left) { return (:: cos((double)_Left)); }
extern "C"    double __cdecl cosh(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type cosh(_Ty _Left) { return (:: cosh((double)_Left)); }
extern "C"    double __cdecl exp(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type exp(_Ty _Left) { return (:: exp((double)_Left)); }

extern "C"    double __cdecl fabs(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type fabs(_Ty _Left) { return (:: fabs((double)_Left)); }

extern "C"   __declspec(dllimport) double __cdecl floor(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type floor(_Ty _Left) { return (:: floor((double)_Left)); }
extern "C"    double __cdecl fmod(   double,   double); template<class _Ty1, class _Ty2> inline typename ::std:: enable_if< ::std:: is_arithmetic<_Ty1>::value && ::std:: is_arithmetic<_Ty2>::value, typename ::std:: _Common_float_type<_Ty1, _Ty2>::type>::type fmod(_Ty1 _Left, _Ty2 _Right) { typedef typename ::std:: _Common_float_type<_Ty1, _Ty2>::type type; return (:: fmod((type)_Left, (type)_Right)); }
extern "C"   __declspec(dllimport) double __cdecl frexp(  double,   int *); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type frexp(_Ty _Left,   int * _Arg2) { return (:: frexp((double)_Left, _Arg2)); }
extern "C"   __declspec(dllimport) double __cdecl ldexp(  double,   int); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type ldexp(_Ty _Left,   int _Arg2) { return (:: ldexp((double)_Left, _Arg2)); }
extern "C"    double __cdecl log(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type log(_Ty _Left) { return (:: log((double)_Left)); }
extern "C"    double __cdecl log10(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type log10(_Ty _Left) { return (:: log10((double)_Left)); }


extern "C"    double __cdecl sin(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type sin(_Ty _Left) { return (:: sin((double)_Left)); }
extern "C"    double __cdecl sinh(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type sinh(_Ty _Left) { return (:: sinh((double)_Left)); }
extern "C"    double __cdecl sqrt(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type sqrt(_Ty _Left) { return (:: sqrt((double)_Left)); }
extern "C"    double __cdecl tan(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type tan(_Ty _Left) { return (:: tan((double)_Left)); }
extern "C"    double __cdecl tanh(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type tanh(_Ty _Left) { return (:: tanh((double)_Left)); }

		









		

inline float _Fma(float _Left, float _Middle, float _Right)
	{	
	return (:: fmaf(_Left, _Middle, _Right));
	}

inline double _Fma(double _Left, double _Middle, double _Right)
	{	
	return (:: fma(_Left, _Middle, _Right));
	}

inline long double _Fma(long double _Left, long double _Middle,
	long double _Right)
	{	
	return (:: fmal(_Left, _Middle, _Right));
	}

template<class _Ty1,
	class _Ty2,
	class _Ty3> inline
	typename ::std:: _Common_float_type<_Ty1,
		typename ::std:: _Common_float_type<_Ty2, _Ty3>::type>::type
	fma(_Ty1 _Left, _Ty2 _Middle, _Ty3 _Right)
	{	
	typedef typename ::std:: _Common_float_type<_Ty1,
		typename ::std:: _Common_float_type<_Ty2, _Ty3>::type>::type type;
	return (_Fma((type)_Left, (type)_Middle, (type)_Right));
	}

		

inline float _Remquo(float _Left, float _Right, int *_Pquo)
	{	
	return (:: remquof(_Left, _Right, _Pquo));
	}

inline double _Remquo(double _Left, double _Right, int *_Pquo)
	{	
	return (:: remquo(_Left, _Right, _Pquo));
	}

inline long double _Remquo(long double _Left, long double _Right, int *_Pquo)
	{	
	return (:: remquol(_Left, _Right, _Pquo));
	}

template<class _Ty1,
	class _Ty2> inline
	typename ::std:: _Common_float_type<_Ty1, _Ty2>::type
	remquo(_Ty1 _Left, _Ty2 _Right, int *_Pquo)
	{	
	typedef typename ::std:: _Common_float_type<_Ty1, _Ty2>::type type;
	return (_Remquo((type)_Left, (type)_Right, _Pquo));
	}

extern "C"   __declspec(dllimport) double __cdecl acosh(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type acosh(_Ty _Left) { return (:: acosh((double)_Left)); }
extern "C"   __declspec(dllimport) double __cdecl asinh(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type asinh(_Ty _Left) { return (:: asinh((double)_Left)); }
extern "C"   __declspec(dllimport) double __cdecl atanh(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type atanh(_Ty _Left) { return (:: atanh((double)_Left)); }
extern "C"   __declspec(dllimport) double __cdecl cbrt(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type cbrt(_Ty _Left) { return (:: cbrt((double)_Left)); }
extern "C"   __declspec(dllimport) double __cdecl copysign(   double,   double); template<class _Ty1, class _Ty2> inline typename ::std:: enable_if< ::std:: is_arithmetic<_Ty1>::value && ::std:: is_arithmetic<_Ty2>::value, typename ::std:: _Common_float_type<_Ty1, _Ty2>::type>::type copysign(_Ty1 _Left, _Ty2 _Right) { typedef typename ::std:: _Common_float_type<_Ty1, _Ty2>::type type; return (:: copysign((type)_Left, (type)_Right)); }
extern "C"   __declspec(dllimport) double __cdecl erf(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type erf(_Ty _Left) { return (:: erf((double)_Left)); }
extern "C"   __declspec(dllimport) double __cdecl erfc(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type erfc(_Ty _Left) { return (:: erfc((double)_Left)); }
extern "C"   __declspec(dllimport) double __cdecl expm1(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type expm1(_Ty _Left) { return (:: expm1((double)_Left)); }
extern "C"   __declspec(dllimport) double __cdecl exp2(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type exp2(_Ty _Left) { return (:: exp2((double)_Left)); }
extern "C"   __declspec(dllimport) double __cdecl fdim(   double,   double); template<class _Ty1, class _Ty2> inline typename ::std:: enable_if< ::std:: is_arithmetic<_Ty1>::value && ::std:: is_arithmetic<_Ty2>::value, typename ::std:: _Common_float_type<_Ty1, _Ty2>::type>::type fdim(_Ty1 _Left, _Ty2 _Right) { typedef typename ::std:: _Common_float_type<_Ty1, _Ty2>::type type; return (:: fdim((type)_Left, (type)_Right)); }

extern "C"   __declspec(dllimport) double __cdecl fmax(   double,   double); template<class _Ty1, class _Ty2> inline typename ::std:: enable_if< ::std:: is_arithmetic<_Ty1>::value && ::std:: is_arithmetic<_Ty2>::value, typename ::std:: _Common_float_type<_Ty1, _Ty2>::type>::type fmax(_Ty1 _Left, _Ty2 _Right) { typedef typename ::std:: _Common_float_type<_Ty1, _Ty2>::type type; return (:: fmax((type)_Left, (type)_Right)); }
extern "C"   __declspec(dllimport) double __cdecl fmin(   double,   double); template<class _Ty1, class _Ty2> inline typename ::std:: enable_if< ::std:: is_arithmetic<_Ty1>::value && ::std:: is_arithmetic<_Ty2>::value, typename ::std:: _Common_float_type<_Ty1, _Ty2>::type>::type fmin(_Ty1 _Left, _Ty2 _Right) { typedef typename ::std:: _Common_float_type<_Ty1, _Ty2>::type type; return (:: fmin((type)_Left, (type)_Right)); }
extern "C"   __declspec(dllimport) double __cdecl hypot(   double,   double); template<class _Ty1, class _Ty2> inline typename ::std:: enable_if< ::std:: is_arithmetic<_Ty1>::value && ::std:: is_arithmetic<_Ty2>::value, typename ::std:: _Common_float_type<_Ty1, _Ty2>::type>::type hypot(_Ty1 _Left, _Ty2 _Right) { typedef typename ::std:: _Common_float_type<_Ty1, _Ty2>::type type; return (:: hypot((type)_Left, (type)_Right)); }
extern "C"   __declspec(dllimport) int __cdecl ilogb(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, int>::type ilogb(_Ty _Left) { return (:: ilogb((double)_Left)); }
extern "C"   __declspec(dllimport) double __cdecl lgamma(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type lgamma(_Ty _Left) { return (:: lgamma((double)_Left)); }
extern "C"   __declspec(dllimport) long long __cdecl llrint(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, long long>::type llrint(_Ty _Left) { return (:: llrint((double)_Left)); }
extern "C"   __declspec(dllimport) long long __cdecl llround(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, long long>::type llround(_Ty _Left) { return (:: llround((double)_Left)); }
extern "C"   __declspec(dllimport) double __cdecl log1p(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type log1p(_Ty _Left) { return (:: log1p((double)_Left)); }
extern "C"   __declspec(dllimport) double __cdecl log2(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type log2(_Ty _Left) { return (:: log2((double)_Left)); }
extern "C"   __declspec(dllimport) double __cdecl logb(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type logb(_Ty _Left) { return (:: logb((double)_Left)); }
extern "C"   __declspec(dllimport) long __cdecl lrint(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, long>::type lrint(_Ty _Left) { return (:: lrint((double)_Left)); }
extern "C"   __declspec(dllimport) long __cdecl lround(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, long>::type lround(_Ty _Left) { return (:: lround((double)_Left)); }
extern "C"   __declspec(dllimport) double __cdecl nearbyint(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type nearbyint(_Ty _Left) { return (:: nearbyint((double)_Left)); }
extern "C"   __declspec(dllimport) double __cdecl nextafter(   double,   double); template<class _Ty1, class _Ty2> inline typename ::std:: enable_if< ::std:: is_arithmetic<_Ty1>::value && ::std:: is_arithmetic<_Ty2>::value, typename ::std:: _Common_float_type<_Ty1, _Ty2>::type>::type nextafter(_Ty1 _Left, _Ty2 _Right) { typedef typename ::std:: _Common_float_type<_Ty1, _Ty2>::type type; return (:: nextafter((type)_Left, (type)_Right)); }
extern "C"   __declspec(dllimport) double __cdecl nexttoward(  double,   long double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type nexttoward(_Ty _Left,   long double _Arg2) { return (:: nexttoward((double)_Left, _Arg2)); }
extern "C"   __declspec(dllimport) double __cdecl remainder(   double,   double); template<class _Ty1, class _Ty2> inline typename ::std:: enable_if< ::std:: is_arithmetic<_Ty1>::value && ::std:: is_arithmetic<_Ty2>::value, typename ::std:: _Common_float_type<_Ty1, _Ty2>::type>::type remainder(_Ty1 _Left, _Ty2 _Right) { typedef typename ::std:: _Common_float_type<_Ty1, _Ty2>::type type; return (:: remainder((type)_Left, (type)_Right)); }

extern "C"   __declspec(dllimport) double __cdecl rint(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type rint(_Ty _Left) { return (:: rint((double)_Left)); }
extern "C"   __declspec(dllimport) double __cdecl round(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type round(_Ty _Left) { return (:: round((double)_Left)); }
extern "C"   __declspec(dllimport) double __cdecl scalbln(  double,   long); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type scalbln(_Ty _Left,   long _Arg2) { return (:: scalbln((double)_Left, _Arg2)); }
extern "C"   __declspec(dllimport) double __cdecl scalbn(  double,   int); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type scalbn(_Ty _Left,   int _Arg2) { return (:: scalbn((double)_Left, _Arg2)); }
extern "C"   __declspec(dllimport) double __cdecl tgamma(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type tgamma(_Ty _Left) { return (:: tgamma((double)_Left)); }
extern "C"   __declspec(dllimport) double __cdecl trunc(  double); template<class _Ty> inline typename ::std:: enable_if< ::std:: is_integral<_Ty>::value, double>::type trunc(_Ty _Left) { return (:: trunc((double)_Left)); }

 
 #pragma warning(pop)
 #pragma pack(pop)


 








 
namespace std {
using :: abs; using :: acos; using :: asin;
using :: atan; using :: atan2; using :: ceil;
using :: cos; using :: cosh; using :: exp;
using :: fabs; using :: floor; using :: fmod;
using :: frexp; using :: ldexp; using :: log;
using :: log10; using :: modf; using :: pow;
using :: sin; using :: sinh; using :: sqrt;
using :: tan; using :: tanh;

using :: acosf; using :: asinf;
using :: atanf; using :: atan2f; using :: ceilf;
using :: cosf; using :: coshf; using :: expf;
using :: fabsf; using :: floorf; using :: fmodf;
using :: frexpf; using :: ldexpf; using :: logf;
using :: log10f; using :: modff; using :: powf;
using :: sinf; using :: sinhf; using :: sqrtf;
using :: tanf; using :: tanhf;

using :: acosl; using :: asinl;
using :: atanl; using :: atan2l; using :: ceill;
using :: cosl; using :: coshl; using :: expl;
using :: fabsl; using :: floorl; using :: fmodl;
using :: frexpl; using :: ldexpl; using :: logl;
using :: log10l; using :: modfl; using :: powl;
using :: sinl; using :: sinhl; using :: sqrtl;
using :: tanl; using :: tanhl;

using :: float_t; using :: double_t;

using :: acosh; using :: asinh; using :: atanh;
using :: cbrt; using :: erf; using :: erfc;
using :: expm1; using :: exp2;
using :: hypot; using :: ilogb; using :: lgamma;
using :: log1p; using :: log2; using :: logb;
using :: llrint; using :: lrint; using :: nearbyint;
using :: rint; using :: llround; using :: lround;
using :: fdim; using :: fma; using :: fmax; using :: fmin;
using :: round; using :: trunc;
using :: remainder; using :: remquo;
using :: copysign; using :: nan; using :: nextafter;
using :: scalbn; using :: scalbln;
using :: nexttoward; using :: tgamma;

using :: acoshf; using :: asinhf; using :: atanhf;
using :: cbrtf; using :: erff; using :: erfcf;
using :: expm1f; using :: exp2f;
using :: hypotf; using :: ilogbf; using :: lgammaf;
using :: log1pf; using :: log2f; using :: logbf;
using :: llrintf; using :: lrintf; using :: nearbyintf;
using :: rintf; using :: llroundf; using :: lroundf;
using :: fdimf; using :: fmaf; using :: fmaxf; using :: fminf;
using :: roundf; using :: truncf;
using :: remainderf; using :: remquof;
using :: copysignf; using :: nanf;
using :: nextafterf; using :: scalbnf; using :: scalblnf;
using :: nexttowardf; using :: tgammaf;

using :: acoshl; using :: asinhl; using :: atanhl;
using :: cbrtl; using :: erfl; using :: erfcl;
using :: expm1l; using :: exp2l;
using :: hypotl; using :: ilogbl; using :: lgammal;
using :: log1pl; using :: log2l; using :: logbl;
using :: llrintl; using :: lrintl; using :: nearbyintl;
using :: rintl; using :: llroundl; using :: lroundl;
using :: fdiml; using :: fmal; using :: fmaxl; using :: fminl;
using :: roundl; using :: truncl;
using :: remainderl; using :: remquol;
using :: copysignl; using :: nanl;
using :: nextafterl; using :: scalbnl; using :: scalblnl;
using :: nexttowardl; using :: tgammal;

using :: fpclassify; using :: signbit;
using :: isfinite; using :: isinf;
using :: isnan; using :: isnormal;
using :: isgreater; using :: isgreaterequal;
using :: isless; using :: islessequal;
using :: islessgreater; using :: isunordered;
}
 










#pragma once










 









#pragma once











#pragma once










#pragma once




__pragma(pack(push, 8)) extern "C" {



__declspec(dllimport) extern int* __cdecl _errno(void);


__declspec(dllimport) errno_t __cdecl _set_errno(  int _Value);
__declspec(dllimport) errno_t __cdecl _get_errno(  int* _Value);



__declspec(dllimport) unsigned long* __cdecl __doserrno(void);


__declspec(dllimport) errno_t __cdecl _set_doserrno(  unsigned long _Value);
__declspec(dllimport) errno_t __cdecl _get_doserrno(  unsigned long * _Value);










































    
    
    
    
    







    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    




} __pragma(pack(pop))









#pragma once





































































































































































































































































































































__pragma(pack(push, 8)) extern "C" {



 
__declspec(dllimport) void const* __cdecl memchr(
      void const* _Buf,
                                 int         _Val,
                                 size_t      _MaxCount
    );

 
int __cdecl memcmp(
      void const* _Buf1,
      void const* _Buf2,
                         size_t      _Size
    );


 

void* __cdecl memcpy(
      void* _Dst,
            void const* _Src,
                               size_t      _Size
    );


__declspec(dllimport) void* __cdecl memmove(
      void*       _Dst,
            void const* _Src,
                                   size_t      _Size
    );

 

void* __cdecl memset(
      void*  _Dst,
                               int    _Val,
                               size_t _Size
    );

 
__declspec(dllimport) char const* __cdecl strchr(
      char const* _Str,
        int         _Val
    );

 
__declspec(dllimport) char const* __cdecl strrchr(
      char const* _Str,
        int         _Ch
    );

   
__declspec(dllimport) char const* __cdecl strstr(
      char const* _Str,
      char const* _SubStr
    );

 

__declspec(dllimport) wchar_t const* __cdecl wcschr(
      wchar_t const* _Str,
        wchar_t        _Ch
    );

 
__declspec(dllimport) wchar_t const* __cdecl wcsrchr(
      wchar_t const* _Str,
        wchar_t        _Ch
    );

   

__declspec(dllimport) wchar_t const* __cdecl wcsstr(
      wchar_t const* _Str,
      wchar_t const* _SubStr
    );



} __pragma(pack(pop))




__pragma(pack(push, 8)) extern "C" {


    















     
    
    static __inline errno_t __cdecl memcpy_s(
          void*       const _Destination,
                                                              rsize_t     const _DestinationSize,
                                 void const* const _Source,
                                                              rsize_t     const _SourceSize
        )
    {
        if (_SourceSize == 0)
        {
            return 0;
        }

        { int _Expr_val=!!(_Destination != 0); if (!(_Expr_val)) { (*_errno()) = 22; _invalid_parameter_noinfo(); return 22; } };
        if (_Source == 0 || _DestinationSize < _SourceSize)
        {
            memset(_Destination, 0, _DestinationSize);

            { int _Expr_val=!!(_Source != 0); if (!(_Expr_val)) { (*_errno()) = 22; _invalid_parameter_noinfo(); return 22; } };
            { int _Expr_val=!!(_DestinationSize >= _SourceSize); if (!(_Expr_val)) { (*_errno()) = 34; _invalid_parameter_noinfo(); return 34; } };

            
            return 22;
        }

        memcpy(_Destination, _Source, _SourceSize);
        return 0;
    }

    
    static __inline errno_t __cdecl memmove_s(
          void*       const _Destination,
                                                              rsize_t     const _DestinationSize,
                                 void const* const _Source,
                                                              rsize_t     const _SourceSize
        )
    {
        if (_SourceSize == 0)
        {
            return 0;
        }

        { int _Expr_val=!!(_Destination != 0); if (!(_Expr_val)) { (*_errno()) = 22; _invalid_parameter_noinfo(); return 22; } };
        { int _Expr_val=!!(_Source != 0); if (!(_Expr_val)) { (*_errno()) = 22; _invalid_parameter_noinfo(); return 22; } };
        { int _Expr_val=!!(_DestinationSize >= _SourceSize); if (!(_Expr_val)) { (*_errno()) = 34; _invalid_parameter_noinfo(); return 34; } };

        memmove(_Destination, _Source, _SourceSize);
        return 0;
    }





} __pragma(pack(pop))










#pragma once










#pragma once





__pragma(pack(push, 8)) extern "C" {






    





    


































    


        #pragma detect_mismatch("_CRT_STDIO_ISO_WIDE_SPECIFIERS", "0")
    




   
__declspec(noinline) __inline unsigned __int64* __cdecl __local_stdio_printf_options(void)
{
    static unsigned __int64 _OptionsStorage;
    return &_OptionsStorage;
}



   
__declspec(noinline) __inline unsigned __int64* __cdecl __local_stdio_scanf_options(void)
{
    static unsigned __int64 _OptionsStorage;
    return &_OptionsStorage;
}



















} __pragma(pack(pop))




__pragma(pack(push, 8)) extern "C" {









    
     
    __declspec(dllimport) errno_t __cdecl _cgetws_s(
          wchar_t* _Buffer,
                                               size_t   _BufferCount,
                                              size_t*  _SizeRead
        );

    extern "C++" { template <size_t _Size> inline   errno_t __cdecl _cgetws_s(  wchar_t (&_Buffer)[_Size],   size_t* _SizeRead) throw() { return _cgetws_s(_Buffer, _Size, _SizeRead); } }

    
    __declspec(dllimport) int __cdecl _cputws(
          wchar_t const* _Buffer
        );

          __declspec(dllimport) wint_t __cdecl _getwch  (void);
          __declspec(dllimport) wint_t __cdecl _getwche (void);
     __declspec(dllimport) wint_t __cdecl _putwch  (  wchar_t _Character);
     __declspec(dllimport) wint_t __cdecl _ungetwch(  wint_t  _Character);

          __declspec(dllimport) wint_t __cdecl _getwch_nolock  (void);
          __declspec(dllimport) wint_t __cdecl _getwche_nolock (void);
     __declspec(dllimport) wint_t __cdecl _putwch_nolock  (  wchar_t _Character);
     __declspec(dllimport) wint_t __cdecl _ungetwch_nolock(  wint_t  _Character);



    
    
    
    
    
    
    __declspec(dllimport) int __cdecl __conio_common_vcwprintf(
                                             unsigned __int64 _Options,
            wchar_t const*   _Format,
                                         _locale_t        _Locale,
                                                va_list          _ArgList
        );

    
    __declspec(dllimport) int __cdecl __conio_common_vcwprintf_s(
                                             unsigned __int64 _Options,
            wchar_t const*   _Format,
                                         _locale_t        _Locale,
                                                va_list          _ArgList
        );

    
    __declspec(dllimport) int __cdecl __conio_common_vcwprintf_p(
                                             unsigned __int64 _Options,
            wchar_t const*   _Format,
                                         _locale_t        _Locale,
                                                va_list          _ArgList
        );

    
    __inline int __cdecl _vcwprintf_l(
            wchar_t const* const _Format,
                                         _locale_t      const _Locale,
                                                va_list              _ArgList
        )



    {
        return __conio_common_vcwprintf((*__local_stdio_printf_options()), _Format, _Locale, _ArgList);
    }


    
    __inline int __cdecl _vcwprintf(
            wchar_t const* const _Format,
                                      va_list              _ArgList
        )



    {
        return _vcwprintf_l(_Format, 0, _ArgList);
    }


    
    __inline int __cdecl _vcwprintf_s_l(
            wchar_t const* const _Format,
                                         _locale_t      const _Locale,
                                                va_list              _ArgList
        )



    {
        return __conio_common_vcwprintf_s((*__local_stdio_printf_options()), _Format, _Locale, _ArgList);
    }


    
    __inline int __cdecl _vcwprintf_s(
            wchar_t const* const _Format,
                                      va_list              _ArgList
        )



    {
        return _vcwprintf_s_l(_Format, 0, _ArgList);
    }


    
    __inline int __cdecl _vcwprintf_p_l(
            wchar_t const* const _Format,
                                         _locale_t      const _Locale,
                                                va_list              _ArgList
        )



    {
        return __conio_common_vcwprintf_p((*__local_stdio_printf_options()), _Format, _Locale, _ArgList);
    }


    
    __inline int __cdecl _vcwprintf_p(
            const wchar_t* const _Format,
                                      va_list              _ArgList
        )



    {
        return _vcwprintf_p_l(_Format, 0, _ArgList);
    }


    
    __inline int __cdecl _cwprintf_l(
            wchar_t const* const _Format,
                                         _locale_t      const _Locale,
        ...)



    {
        int _Result;
        va_list _ArgList;
        ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
        _Result = _vcwprintf_l(_Format, _Locale, _ArgList);
        ((void)(_ArgList = (va_list)0));
        return _Result;
    }


    
    __inline int __cdecl _cwprintf(
            wchar_t const* const _Format,
        ...)



    {
        int _Result;
        va_list _ArgList;
        ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
        _Result = _vcwprintf_l(_Format, 0, _ArgList);
        ((void)(_ArgList = (va_list)0));
        return _Result;
    }


    
    __inline int __cdecl _cwprintf_s_l(
            wchar_t const* const _Format,
                                         _locale_t      const _Locale,
        ...)



    {
        int _Result;
        va_list _ArgList;
        ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
        _Result = _vcwprintf_s_l(_Format, _Locale, _ArgList);
        ((void)(_ArgList = (va_list)0));
        return _Result;
    }


    
    __inline int __cdecl _cwprintf_s(
            wchar_t const* const _Format,
        ...)



    {
        int _Result;
        va_list _ArgList;
        ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
        _Result = _vcwprintf_s_l(_Format, 0, _ArgList);
        ((void)(_ArgList = (va_list)0));
        return _Result;
    }


    
    __inline int __cdecl _cwprintf_p_l(
            wchar_t const* const _Format,
                                         _locale_t      const _Locale,
        ...)



    {
        int _Result;
        va_list _ArgList;
        ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
        _Result = _vcwprintf_p_l(_Format, _Locale, _ArgList);
        ((void)(_ArgList = (va_list)0));
        return _Result;
    }


    
    __inline int __cdecl _cwprintf_p(
            wchar_t const* const _Format,
        ...)



    {
        int _Result;
        va_list _ArgList;
        ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
        _Result = _vcwprintf_p_l(_Format, 0, _ArgList);
        ((void)(_ArgList = (va_list)0));
        return _Result;
    }




    
    
    
    
    
    
    __declspec(dllimport) int __cdecl __conio_common_vcwscanf(
                                            unsigned __int64 _Options,
            wchar_t const*   _Format,
                                        _locale_t        _Locale,
                                               va_list          _ArgList
        );

     __declspec(deprecated("This function or variable may be unsafe. Consider using " "_vcwscanf_s_l" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
    __inline int __cdecl _vcwscanf_l(
            wchar_t const* const _Format,
                                        _locale_t      const _Locale,
                                               va_list              _ArgList
        )



    {
        return __conio_common_vcwscanf(
            (*__local_stdio_scanf_options ()),
            _Format, _Locale, _ArgList);
    }


     __declspec(deprecated("This function or variable may be unsafe. Consider using " "_vcwscanf_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
    __inline int __cdecl _vcwscanf(
            wchar_t const* const _Format,
                                               va_list              _ArgList
        )



    {
        #pragma warning(push)
        #pragma warning(disable: 4996) 
        return _vcwscanf_l(_Format, 0, _ArgList);
        #pragma warning(pop)
    }


    
    __inline int __cdecl _vcwscanf_s_l(
            wchar_t const* const _Format,
                                        _locale_t      const _Locale,
                                               va_list              _ArgList
        )



    {
        return __conio_common_vcwscanf(
            (*__local_stdio_scanf_options ()) | (1ULL << 0),
            _Format, _Locale, _ArgList);
    }


    
    __inline int __cdecl _vcwscanf_s(
            wchar_t const* const _Format,
                                               va_list              _ArgList
        )



    {
        return _vcwscanf_s_l(_Format, 0, _ArgList);
    }


     __declspec(deprecated("This function or variable may be unsafe. Consider using " "_cwscanf_s_l" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
    __inline int __cdecl _cwscanf_l(
            wchar_t const* const _Format,
                                        _locale_t      const _Locale,
        ...)



    {
        int _Result;
        va_list _ArgList;
        ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));

        #pragma warning(push)
        #pragma warning(disable: 4996) 
        _Result = _vcwscanf_l(_Format, _Locale, _ArgList);
        #pragma warning(pop)

        ((void)(_ArgList = (va_list)0));
        return _Result;
    }


     __declspec(deprecated("This function or variable may be unsafe. Consider using " "_cwscanf_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
    __inline int __cdecl _cwscanf(
            wchar_t const* const _Format,
        ...)



    {
        int _Result;
        va_list _ArgList;
        ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));

        #pragma warning(push)
        #pragma warning(disable: 4996) 
        _Result = _vcwscanf_l(_Format, 0, _ArgList);
        #pragma warning(pop)

        ((void)(_ArgList = (va_list)0));
        return _Result;
    }


    
    __inline int __cdecl _cwscanf_s_l(
            wchar_t const* const _Format,
                                        _locale_t      const _Locale,
        ...)



    {
        int _Result;
        va_list _ArgList;
        ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
        _Result = _vcwscanf_s_l(_Format, _Locale, _ArgList);
        ((void)(_ArgList = (va_list)0));
        return _Result;
    }


    
    __inline int __cdecl _cwscanf_s(
            wchar_t const* const _Format,
        ...)



    {
        int _Result;
        va_list _ArgList;
        ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
        _Result = _vcwscanf_s_l(_Format, 0, _ArgList);
        ((void)(_ArgList = (va_list)0));
        return _Result;
    }






} __pragma(pack(pop))










#pragma once



__pragma(pack(push, 8)) extern "C" {











    


        
    

    __declspec(dllimport) const unsigned short* __cdecl __pctype_func(void);
    __declspec(dllimport) const wctype_t*       __cdecl __pwctype_func(void);

    



        
        
    






















  __declspec(dllimport) int __cdecl iswalnum  (  wint_t _C);
  __declspec(dllimport) int __cdecl iswalpha  (  wint_t _C);
  __declspec(dllimport) int __cdecl iswascii  (  wint_t _C);
  __declspec(dllimport) int __cdecl iswblank  (  wint_t _C);
  __declspec(dllimport) int __cdecl iswcntrl  (  wint_t _C);


  __declspec(dllimport) int __cdecl iswdigit  (  wint_t _C);

  __declspec(dllimport) int __cdecl iswgraph  (  wint_t _C);
  __declspec(dllimport) int __cdecl iswlower  (  wint_t _C);
  __declspec(dllimport) int __cdecl iswprint  (  wint_t _C);
  __declspec(dllimport) int __cdecl iswpunct  (  wint_t _C);
  __declspec(dllimport) int __cdecl iswspace  (  wint_t _C);
  __declspec(dllimport) int __cdecl iswupper  (  wint_t _C);
  __declspec(dllimport) int __cdecl iswxdigit (  wint_t _C);
  __declspec(dllimport) int __cdecl __iswcsymf(  wint_t _C);
  __declspec(dllimport) int __cdecl __iswcsym (  wint_t _C);

  __declspec(dllimport) int __cdecl _iswalnum_l (  wint_t _C,   _locale_t _Locale);
  __declspec(dllimport) int __cdecl _iswalpha_l (  wint_t _C,   _locale_t _Locale);
  __declspec(dllimport) int __cdecl _iswblank_l (  wint_t _C,   _locale_t _Locale);
  __declspec(dllimport) int __cdecl _iswcntrl_l (  wint_t _C,   _locale_t _Locale);
  __declspec(dllimport) int __cdecl _iswdigit_l (  wint_t _C,   _locale_t _Locale);
  __declspec(dllimport) int __cdecl _iswgraph_l (  wint_t _C,   _locale_t _Locale);
  __declspec(dllimport) int __cdecl _iswlower_l (  wint_t _C,   _locale_t _Locale);
  __declspec(dllimport) int __cdecl _iswprint_l (  wint_t _C,   _locale_t _Locale);
  __declspec(dllimport) int __cdecl _iswpunct_l (  wint_t _C,   _locale_t _Locale);
  __declspec(dllimport) int __cdecl _iswspace_l (  wint_t _C,   _locale_t _Locale);
  __declspec(dllimport) int __cdecl _iswupper_l (  wint_t _C,   _locale_t _Locale);
  __declspec(dllimport) int __cdecl _iswxdigit_l(  wint_t _C,   _locale_t _Locale);
  __declspec(dllimport) int __cdecl _iswcsymf_l (  wint_t _C,   _locale_t _Locale);
  __declspec(dllimport) int __cdecl _iswcsym_l  (  wint_t _C,   _locale_t _Locale);


  __declspec(dllimport) wint_t __cdecl towupper(  wint_t _C);
  __declspec(dllimport) wint_t __cdecl towlower(  wint_t _C);
  __declspec(dllimport) int    __cdecl iswctype(  wint_t _C,   wctype_t _Type);

  __declspec(dllimport) wint_t __cdecl _towupper_l(  wint_t _C,   _locale_t _Locale);
  __declspec(dllimport) wint_t __cdecl _towlower_l(  wint_t _C,   _locale_t _Locale);
  __declspec(dllimport) int    __cdecl _iswctype_l(  wint_t _C,   wctype_t _Type,   _locale_t _Locale);



      __declspec(dllimport) int __cdecl isleadbyte(  int _C);
      __declspec(dllimport) int __cdecl _isleadbyte_l(  int _C,   _locale_t _Locale);

    __declspec(deprecated("This function or variable has been superceded by newer library " "or operating system functionality. Consider using " "iswctype" " " "instead. See online help for details.")) __declspec(dllimport) int __cdecl is_wctype(  wint_t _C,   wctype_t _Type);























































































} __pragma(pack(pop))










#pragma once



__pragma(pack(push, 8)) extern "C" {






 
   
__declspec(dllimport) __declspec(allocator) wchar_t* __cdecl _wgetcwd(
      wchar_t* _DstBuf,
                                  int      _SizeInWords
    );

 
   
__declspec(dllimport) __declspec(allocator) wchar_t* __cdecl _wgetdcwd(
                                  int      _Drive,
      wchar_t* _DstBuf,
                                  int      _SizeInWords
    );






 
__declspec(dllimport) int __cdecl _wchdir(
      wchar_t const* _Path
    );

 
__declspec(dllimport) int __cdecl _wmkdir(
      wchar_t const* _Path
    );

 
__declspec(dllimport) int __cdecl _wrmdir(
      wchar_t const* _Path
    );



} __pragma(pack(pop))










#pragma once











#pragma once












    
    
    
    



__pragma(pack(push, 8)) extern "C" {



#pragma warning(push)
#pragma warning(disable:4820) 












    
    


typedef unsigned long _fsize_t;

struct _wfinddata32_t
{
    unsigned   attrib;
    __time32_t time_create;    
    __time32_t time_access;    
    __time32_t time_write;
    _fsize_t   size;
    wchar_t    name[260];
};

struct _wfinddata32i64_t
{
    unsigned   attrib;
    __time32_t time_create;    
    __time32_t time_access;    
    __time32_t time_write;
    __int64    size;
    wchar_t    name[260];
};

struct _wfinddata64i32_t
{
    unsigned   attrib;
    __time64_t time_create;    
    __time64_t time_access;    
    __time64_t time_write;
    _fsize_t   size;
    wchar_t    name[260];
};

struct _wfinddata64_t
{
    unsigned   attrib;
    __time64_t time_create;    
    __time64_t time_access;    
    __time64_t time_write;
    __int64    size;
    wchar_t    name[260];
};














    
    
    
    


 
__declspec(dllimport) int __cdecl _waccess(
      wchar_t const* _FileName,
        int            _AccessMode
    );


__declspec(dllimport) errno_t __cdecl _waccess_s(
      wchar_t const* _FileName,
        int            _AccessMode
    );

 
__declspec(dllimport) int __cdecl _wchmod(
      wchar_t const* _FileName,
        int            _Mode
    );

  __declspec(deprecated("This function or variable may be unsafe. Consider using " "_wsopen_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__declspec(dllimport) int __cdecl _wcreat(
      wchar_t const* _FileName,
        int            _PermissionMode
    );

 
 
__declspec(dllimport) intptr_t __cdecl _wfindfirst32(
      wchar_t const*         _FileName,
       struct _wfinddata32_t* _FindData
    );

 
 
__declspec(dllimport) int __cdecl _wfindnext32(
       intptr_t               _FindHandle,
      struct _wfinddata32_t* _FindData
    );

__declspec(dllimport) int __cdecl _wunlink(
      wchar_t const* _FileName
    );

 
__declspec(dllimport) int __cdecl _wrename(
      wchar_t const* _OldFileName,
      wchar_t const* _NewFileName
    );

__declspec(dllimport) errno_t __cdecl _wmktemp_s(
      wchar_t* _TemplateName,
                                 size_t   _SizeInWords
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl _wmktemp_s(wchar_t (&_TemplateName)[_Size]) throw() { return _wmktemp_s(_TemplateName, _Size); } }

 
__declspec(deprecated("This function or variable may be unsafe. Consider using " "_wmktemp_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport) wchar_t* __cdecl _wmktemp( wchar_t *_TemplateName);

 
 
__declspec(dllimport) intptr_t __cdecl _wfindfirst32i64(
      wchar_t const*            _FileName,
       struct _wfinddata32i64_t* _FindData
    );

 
 
__declspec(dllimport) intptr_t __cdecl _wfindfirst64i32(
      wchar_t const*            _FileName,
       struct _wfinddata64i32_t* _FindData
    );

 
 
__declspec(dllimport) intptr_t __cdecl _wfindfirst64(
      wchar_t const*         _FileName,
       struct _wfinddata64_t* _FindData
    );

 
 
__declspec(dllimport) int __cdecl _wfindnext32i64(
       intptr_t                  _FindHandle,
      struct _wfinddata32i64_t* _FindData
    );

 
 
__declspec(dllimport) int __cdecl _wfindnext64i32(
       intptr_t                  _FindHandle,
      struct _wfinddata64i32_t* _FindData
    );

 
 
__declspec(dllimport) int __cdecl _wfindnext64(
       intptr_t               _FindHandle,
      struct _wfinddata64_t* _FindData
    );


__declspec(dllimport) errno_t __cdecl _wsopen_s(
       int*           _FileHandle,
      wchar_t const* _FileName,
        int            _OpenFlag,
        int            _ShareFlag,
        int            _PermissionFlag
    );

__declspec(dllimport) errno_t __cdecl _wsopen_dispatch(
      wchar_t const* _FileName,
        int            _OFlag,
        int            _ShFlag,
        int            _PMode,
       int*           _PFileHandle,
        int            _BSecure
    );





    
    extern "C++"   __declspec(deprecated("This function or variable may be unsafe. Consider using " "_wsopen_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
    inline int __cdecl _wopen(
          wchar_t const* _FileName,
            int            _OFlag,
            int            _PMode = 0
        )
    {
        int _FileHandle;
        
        errno_t const _Result = _wsopen_dispatch(_FileName, _OFlag, 0x40, _PMode, &_FileHandle, 0);
        return _Result ? -1 : _FileHandle;
    }

    extern "C++"   __declspec(deprecated("This function or variable may be unsafe. Consider using " "_wsopen_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
    inline int __cdecl _wsopen(
          wchar_t const* _FileName,
            int            _OFlag,
            int            _ShFlag,
            int            _PMode = 0
        )
    {
        int _FileHandle;
        
        errno_t const _Result = _wsopen_dispatch(_FileName, _OFlag, _ShFlag, _PMode, &_FileHandle, 0);
        return _Result ? -1 : _FileHandle;
    }





















#pragma warning(pop)



} __pragma(pack(pop))










#pragma once



__pragma(pack(push, 8)) extern "C" {


    


    __declspec(dllimport) intptr_t __cdecl _wexecl(
          wchar_t const* _FileName,
          wchar_t const* _ArgList,
        ...);

    __declspec(dllimport) intptr_t __cdecl _wexecle(
          wchar_t const* _FileName,
          wchar_t const* _ArgList,
        ...);

    __declspec(dllimport) intptr_t __cdecl _wexeclp(
          wchar_t const* _FileName,
          wchar_t const* _ArgList,
        ...);

    __declspec(dllimport) intptr_t __cdecl _wexeclpe(
          wchar_t const* _FileName,
          wchar_t const* _ArgList,
        ...);

    __declspec(dllimport) intptr_t __cdecl _wexecv(
          wchar_t const*        _FileName,
          wchar_t const* const* _ArgList
        );

    __declspec(dllimport) intptr_t __cdecl _wexecve(
              wchar_t const*        _FileName,
              wchar_t const* const* _ArgList,
          wchar_t const* const* _Env
        );

    __declspec(dllimport) intptr_t __cdecl _wexecvp(
          wchar_t const*        _FileName,
          wchar_t const* const* _ArgList
        );

    __declspec(dllimport) intptr_t __cdecl _wexecvpe(
              wchar_t const*        _FileName,
              wchar_t const* const* _ArgList,
          wchar_t const* const* _Env
        );

    __declspec(dllimport) intptr_t __cdecl _wspawnl(
            int            _Mode,
          wchar_t const* _FileName,
          wchar_t const* _ArgList,
        ...);

    __declspec(dllimport) intptr_t __cdecl _wspawnle(
            int            _Mode,
          wchar_t const* _FileName,
          wchar_t const* _ArgList,
        ...);

    __declspec(dllimport) intptr_t __cdecl _wspawnlp(
            int            _Mode,
          wchar_t const* _FileName,
          wchar_t const* _ArgList,
        ...);

    __declspec(dllimport) intptr_t __cdecl _wspawnlpe(
            int            _Mode,
          wchar_t const* _FileName,
          wchar_t const* _ArgList,
        ...);

    __declspec(dllimport) intptr_t __cdecl _wspawnv(
            int                   _Mode,
          wchar_t const*        _FileName,
          wchar_t const* const* _ArgList
        );

    __declspec(dllimport) intptr_t __cdecl _wspawnve(
                int                   _Mode,
              wchar_t const*        _FileName,
              wchar_t const* const* _ArgList,
          wchar_t const* const* _Env
        );

    __declspec(dllimport) intptr_t __cdecl _wspawnvp(
            int                   _Mode,
          wchar_t const*        _FileName,
          wchar_t const* const* _ArgList
        );

    __declspec(dllimport) intptr_t __cdecl _wspawnvpe(
                int                   _Mode,
              wchar_t const*        _FileName,
              wchar_t const* const* _ArgList,
          wchar_t const* const* _Env
        );

    __declspec(dllimport) int __cdecl _wsystem(
          wchar_t const* _Command
        );





} __pragma(pack(pop))











#pragma once




__pragma(pack(push, 8)) extern "C" {








    
    typedef struct _iobuf
    {
        void* _Placeholder;
    } FILE;


__declspec(dllimport) FILE* __cdecl __acrt_iob_func(unsigned);















__declspec(dllimport) wint_t __cdecl fgetwc(
      FILE* _Stream
    );


__declspec(dllimport) wint_t __cdecl _fgetwchar(void);


__declspec(dllimport) wint_t __cdecl fputwc(
         wchar_t _Character,
      FILE*   _Stream);


__declspec(dllimport) wint_t __cdecl _fputwchar(
      wchar_t _Character
    );

 
__declspec(dllimport) wint_t __cdecl getwc(
      FILE* _Stream
    );

 
__declspec(dllimport) wint_t __cdecl getwchar(void);



 
__declspec(dllimport) wchar_t* __cdecl fgetws(
      wchar_t* _Buffer,
                              int      _BufferCount,
                           FILE*    _Stream
    );


__declspec(dllimport) int __cdecl fputws(
       wchar_t const* _Buffer,
      FILE*          _Stream
    );


 
__declspec(dllimport) wchar_t* __cdecl _getws_s(
      wchar_t* _Buffer,
                              size_t   _BufferCount
    );

extern "C++" { template <size_t _Size> inline   wchar_t* __cdecl _getws_s(  wchar_t (&_Buffer)[_Size]) throw() { return _getws_s(_Buffer, _Size); } }


__declspec(dllimport) wint_t __cdecl putwc(
         wchar_t _Character,
      FILE*   _Stream
    );


__declspec(dllimport) wint_t __cdecl putwchar(
      wchar_t _Character
    );


__declspec(dllimport) int __cdecl _putws(
      wchar_t const* _Buffer
    );


__declspec(dllimport) wint_t __cdecl ungetwc(
         wint_t _Character,
      FILE*  _Stream
    );

 
__declspec(dllimport) FILE * __cdecl _wfdopen(
        int            _FileHandle,
      wchar_t const* _Mode
    );

  __declspec(deprecated("This function or variable may be unsafe. Consider using " "_wfopen_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__declspec(dllimport) FILE* __cdecl _wfopen(
      wchar_t const* _FileName,
      wchar_t const* _Mode
    );


__declspec(dllimport) errno_t __cdecl _wfopen_s(
      FILE**         _Stream,
                         wchar_t const* _FileName,
                         wchar_t const* _Mode
    );

 
__declspec(deprecated("This function or variable may be unsafe. Consider using " "_wfreopen_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__declspec(dllimport) FILE* __cdecl _wfreopen(
       wchar_t const* _FileName,
       wchar_t const* _Mode,
      FILE*          _OldStream
    );


__declspec(dllimport) errno_t __cdecl _wfreopen_s(
      FILE**         _Stream,
                         wchar_t const* _FileName,
                         wchar_t const* _Mode,
                        FILE*          _OldStream
    );

 
__declspec(dllimport) FILE* __cdecl _wfsopen(
      wchar_t const* _FileName,
      wchar_t const* _Mode,
        int            _ShFlag
    );

__declspec(dllimport) void __cdecl _wperror(
      wchar_t const* _ErrorMessage
    );



     
    __declspec(dllimport) FILE* __cdecl _wpopen(
          wchar_t const* _Command,
          wchar_t const* _Mode
        );



__declspec(dllimport) int __cdecl _wremove(
      wchar_t const* _FileName
    );




 
__declspec(dllimport) __declspec(allocator) wchar_t* __cdecl _wtempnam(
      wchar_t const* _Directory,
      wchar_t const* _FilePrefix
    );



 

__declspec(dllimport) errno_t __cdecl _wtmpnam_s(
      wchar_t* _Buffer,
                              size_t   _BufferCount
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl _wtmpnam_s(  wchar_t (&_Buffer)[_Size]) throw() { return _wtmpnam_s(_Buffer, _Size); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "_wtmpnam_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport)  wchar_t* __cdecl _wtmpnam(  wchar_t *_Buffer);









__declspec(dllimport) wint_t __cdecl _fgetwc_nolock(
      FILE* _Stream
    );


__declspec(dllimport) wint_t __cdecl _fputwc_nolock(
         wchar_t _Character, 
      FILE*   _Stream
    );


__declspec(dllimport) wint_t __cdecl _getwc_nolock(
      FILE* _Stream
    );


__declspec(dllimport) wint_t __cdecl _putwc_nolock(
         wchar_t _Character,
      FILE*   _Stream
    );


__declspec(dllimport) wint_t __cdecl _ungetwc_nolock(
         wint_t _Character,
      FILE*  _Stream
    );






















__declspec(dllimport) int __cdecl __stdio_common_vfwprintf(
                                         unsigned __int64 _Options,
                                      FILE*            _Stream,
        wchar_t const*   _Format,
                                     _locale_t        _Locale,
                                            va_list          _ArgList
    );


__declspec(dllimport) int __cdecl __stdio_common_vfwprintf_s(
                                         unsigned __int64 _Options,
                                      FILE*            _Stream,
        wchar_t const*   _Format,
                                     _locale_t        _Locale,
                                            va_list          _ArgList
    );


__declspec(dllimport) int __cdecl __stdio_common_vfwprintf_p(
                                         unsigned __int64 _Options,
                                      FILE*            _Stream,
        wchar_t const*   _Format,
                                     _locale_t        _Locale,
                                            va_list          _ArgList
    );


__inline int __cdecl _vfwprintf_l(
                                      FILE*          const _Stream,
        wchar_t const* const _Format,
                                     _locale_t      const _Locale,
                                            va_list              _ArgList
    )



{
    return __stdio_common_vfwprintf((*__local_stdio_printf_options()), _Stream, _Format, _Locale, _ArgList);
}



__inline int __cdecl vfwprintf(
                            FILE*          const _Stream,
        wchar_t const* const _Format,
                                  va_list              _ArgList
    )



{
    return _vfwprintf_l(_Stream, _Format, 0, _ArgList);
}



__inline int __cdecl _vfwprintf_s_l(
                                      FILE*          const _Stream,
        wchar_t const* const _Format,
                                     _locale_t      const _Locale,
                                            va_list              _ArgList
    )



{
    return __stdio_common_vfwprintf_s((*__local_stdio_printf_options()), _Stream, _Format, _Locale, _ArgList);
}




    
    __inline int __cdecl vfwprintf_s(
                                FILE*          const _Stream,
            wchar_t const* const _Format,
                                      va_list              _ArgList
        )



    {
        return _vfwprintf_s_l(_Stream, _Format, 0, _ArgList);
    }





__inline int __cdecl _vfwprintf_p_l(
                                      FILE*          const _Stream,
        wchar_t const* const _Format,
                                     _locale_t      const _Locale,
                                            va_list              _ArgList
    )



{
    return __stdio_common_vfwprintf_p((*__local_stdio_printf_options()), _Stream, _Format, _Locale, _ArgList);
}



__inline int __cdecl _vfwprintf_p(
                            FILE*          const _Stream,
        wchar_t const* const _Format,
                                  va_list              _ArgList
    )



{
    return _vfwprintf_p_l(_Stream, _Format, 0, _ArgList);
}



__inline int __cdecl _vwprintf_l(
        wchar_t const* const _Format,
                                     _locale_t      const _Locale,
                                            va_list              _ArgList
    )



{
    return _vfwprintf_l((__acrt_iob_func(1)), _Format, _Locale, _ArgList);
}



__inline int __cdecl vwprintf(
        wchar_t const* const _Format,
                                  va_list              _ArgList
    )



{
    return _vfwprintf_l((__acrt_iob_func(1)), _Format, 0, _ArgList);
}



__inline int __cdecl _vwprintf_s_l(
        wchar_t const* const _Format,
                                     _locale_t      const _Locale,
                                            va_list              _ArgList
    )



{
    return _vfwprintf_s_l((__acrt_iob_func(1)), _Format, _Locale, _ArgList);
}




    
    __inline int __cdecl vwprintf_s(
            wchar_t const* const _Format,
                                      va_list              _ArgList
        )



    {
        return _vfwprintf_s_l((__acrt_iob_func(1)), _Format, 0, _ArgList);
    }





__inline int __cdecl _vwprintf_p_l(
        wchar_t const* const _Format,
                                     _locale_t      const _Locale,
                                            va_list              _ArgList
    )



{
    return _vfwprintf_p_l((__acrt_iob_func(1)), _Format, _Locale, _ArgList);
}



__inline int __cdecl _vwprintf_p(
        wchar_t const* const _Format,
                                  va_list              _ArgList
    )



{
    return _vfwprintf_p_l((__acrt_iob_func(1)), _Format, 0, _ArgList);
}



__inline int __cdecl _fwprintf_l(
                                      FILE*          const _Stream,
        wchar_t const* const _Format,
                                     _locale_t      const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vfwprintf_l(_Stream, _Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}



__inline int __cdecl fwprintf(
                            FILE*          const _Stream,
        wchar_t const* const _Format,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
    _Result = _vfwprintf_l(_Stream, _Format, 0, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}



__inline int __cdecl _fwprintf_s_l(
                                      FILE*          const _Stream,
        wchar_t const* const _Format,
                                     _locale_t      const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vfwprintf_s_l(_Stream, _Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}




    
    __inline int __cdecl fwprintf_s(
                                FILE*          const _Stream,
            wchar_t const* const _Format,
        ...)



    {
        int _Result;
        va_list _ArgList;
        ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
        _Result = _vfwprintf_s_l(_Stream, _Format, 0, _ArgList);
        ((void)(_ArgList = (va_list)0));
        return _Result;
    }





__inline int __cdecl _fwprintf_p_l(
                                      FILE*          const _Stream,
        wchar_t const* const _Format,
                                     _locale_t      const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vfwprintf_p_l(_Stream, _Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}



__inline int __cdecl _fwprintf_p(
                            FILE*          const _Stream,
        wchar_t const* const _Format,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
    _Result = _vfwprintf_p_l(_Stream, _Format, 0, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}



__inline int __cdecl _wprintf_l(
        wchar_t const* const _Format,
                                     _locale_t      const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vfwprintf_l((__acrt_iob_func(1)), _Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}



__inline int __cdecl wprintf(
        wchar_t const* const _Format,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
    _Result = _vfwprintf_l((__acrt_iob_func(1)), _Format, 0, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}



__inline int __cdecl _wprintf_s_l(
        wchar_t const* const _Format,
                                     _locale_t      const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vfwprintf_s_l((__acrt_iob_func(1)), _Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}




    
    __inline int __cdecl wprintf_s(
            wchar_t const* const _Format,
        ...)



    {
        int _Result;
        va_list _ArgList;
        ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
        _Result = _vfwprintf_s_l((__acrt_iob_func(1)), _Format, 0, _ArgList);
        ((void)(_ArgList = (va_list)0));
        return _Result;
    }





__inline int __cdecl _wprintf_p_l(
        wchar_t const* const _Format,
                                     _locale_t      const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vfwprintf_p_l((__acrt_iob_func(1)), _Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}



__inline int __cdecl _wprintf_p(
        wchar_t const* const _Format,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
    _Result = _vfwprintf_p_l((__acrt_iob_func(1)), _Format, 0, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}









__declspec(dllimport) int __cdecl __stdio_common_vfwscanf(
                                        unsigned __int64 _Options,
                                     FILE*            _Stream,
        wchar_t const*   _Format,
                                    _locale_t        _Locale,
                                           va_list          _ArgList
    );


__inline int __cdecl _vfwscanf_l(
      FILE*                                const _Stream,
        wchar_t const* const _Format,
                           _locale_t      const _Locale,
                                  va_list              _ArgList
    )



{
    return __stdio_common_vfwscanf(
        (*__local_stdio_scanf_options ()),
        _Stream, _Format, _Locale, _ArgList);
}



__inline int __cdecl vfwscanf(
      FILE*                                const _Stream,
        wchar_t const* const _Format,
                                  va_list              _ArgList
    )



{
    return _vfwscanf_l(_Stream, _Format, 0, _ArgList);
}



__inline int __cdecl _vfwscanf_s_l(
                            FILE*          const _Stream,
        wchar_t const* const _Format,
                           _locale_t      const _Locale,
                                  va_list              _ArgList
    )



{
    return __stdio_common_vfwscanf(
        (*__local_stdio_scanf_options ()) | (1ULL << 0),
        _Stream, _Format, _Locale, _ArgList);
}




    
    __inline int __cdecl vfwscanf_s(
                                FILE*          const _Stream,
            wchar_t const* const _Format,
                                      va_list              _ArgList
        )



    {
        return _vfwscanf_s_l(_Stream, _Format, 0, _ArgList);
    }




__inline int __cdecl _vwscanf_l(
        wchar_t const* const _Format,
                           _locale_t      const _Locale,
                                  va_list              _ArgList
    )



{
    return _vfwscanf_l((__acrt_iob_func(0)), _Format, _Locale, _ArgList);
}



__inline int __cdecl vwscanf(
        wchar_t const* const _Format,
                                  va_list              _ArgList
    )



{
    return _vfwscanf_l((__acrt_iob_func(0)), _Format, 0, _ArgList);
}



__inline int __cdecl _vwscanf_s_l(
        wchar_t const* const _Format,
                           _locale_t      const _Locale,
                                  va_list              _ArgList
    )



{
    return _vfwscanf_s_l((__acrt_iob_func(0)), _Format, _Locale, _ArgList);
}




    
    __inline int __cdecl vwscanf_s(
            wchar_t const* const _Format,
                                      va_list              _ArgList
        )



    {
        return _vfwscanf_s_l((__acrt_iob_func(0)), _Format, 0, _ArgList);
    }




 __declspec(deprecated("This function or variable may be unsafe. Consider using " "_fwscanf_s_l" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__inline int __cdecl _fwscanf_l(
                                     FILE*          const _Stream,
        wchar_t const* const _Format,
                                    _locale_t      const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vfwscanf_l(_Stream, _Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}


  __declspec(deprecated("This function or variable may be unsafe. Consider using " "fwscanf_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__inline int __cdecl fwscanf(
                           FILE*          const _Stream,
        wchar_t const* const _Format,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
    _Result = _vfwscanf_l(_Stream, _Format, 0, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}



__inline int __cdecl _fwscanf_s_l(
                                       FILE*          const _Stream,
        wchar_t const* const _Format,
                                      _locale_t      const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vfwscanf_s_l(_Stream, _Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}




    
    __inline int __cdecl fwscanf_s(
                                 FILE*          const _Stream,
            wchar_t const* const _Format,
        ...)



    {
        int _Result;
        va_list _ArgList;
        ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
        _Result = _vfwscanf_s_l(_Stream, _Format, 0, _ArgList);
        ((void)(_ArgList = (va_list)0));
        return _Result;
    }




 __declspec(deprecated("This function or variable may be unsafe. Consider using " "_wscanf_s_l" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__inline int __cdecl _wscanf_l(
        wchar_t const* const _Format,
                                    _locale_t      const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vfwscanf_l((__acrt_iob_func(0)), _Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}


  __declspec(deprecated("This function or variable may be unsafe. Consider using " "wscanf_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__inline int __cdecl wscanf(
        wchar_t const* const _Format,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
    _Result = _vfwscanf_l((__acrt_iob_func(0)), _Format, 0, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}



__inline int __cdecl _wscanf_s_l(
        wchar_t const* const _Format,
                                      _locale_t      const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vfwscanf_s_l((__acrt_iob_func(0)), _Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}




    
    __inline int __cdecl wscanf_s(
            wchar_t const* const _Format,
        ...)



    {
        int _Result;
        va_list _ArgList;
        ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
        _Result = _vfwscanf_s_l((__acrt_iob_func(0)), _Format, 0, _ArgList);
        ((void)(_ArgList = (va_list)0));
        return _Result;
    }












    







 

__declspec(dllimport) int __cdecl __stdio_common_vswprintf(
                                         unsigned __int64 _Options,
                 wchar_t*         _Buffer,
                                         size_t           _BufferCount,
        wchar_t const*   _Format,
                                     _locale_t        _Locale,
                                            va_list          _ArgList
    );

 

__declspec(dllimport) int __cdecl __stdio_common_vswprintf_s(
                                         unsigned __int64 _Options,
                 wchar_t*         _Buffer,
                                         size_t           _BufferCount,
        wchar_t const*   _Format,
                                     _locale_t        _Locale,
                                            va_list          _ArgList
    );

 

__declspec(dllimport) int __cdecl __stdio_common_vsnwprintf_s(
                                         unsigned __int64 _Options,
                 wchar_t*         _Buffer,
                                         size_t           _BufferCount,
                                         size_t           _MaxCount,
        wchar_t const*   _Format,
                                     _locale_t        _Locale,
                                            va_list          _ArgList
    );

 

__declspec(dllimport) int __cdecl __stdio_common_vswprintf_p(
                                         unsigned __int64 _Options,
                 wchar_t*         _Buffer,
                                         size_t           _BufferCount,
        wchar_t const*   _Format,
                                     _locale_t        _Locale,
                                            va_list          _ArgList
    );

 
 __declspec(deprecated("This function or variable may be unsafe. Consider using " "_vsnwprintf_s_l" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__inline int __cdecl _vsnwprintf_l(
           wchar_t*       const _Buffer,
                                             size_t         const _BufferCount,
            wchar_t const* const _Format,
                                         _locale_t      const _Locale,
                                                va_list              _ArgList
    )



{
    int const _Result = __stdio_common_vswprintf(
        (*__local_stdio_printf_options()) | (1ULL << 0),
        _Buffer, _BufferCount, _Format, _Locale, _ArgList);

    return _Result < 0 ? -1 : _Result;
}


#pragma warning(push)
#pragma warning(disable: 4793)

 

__inline int __cdecl _vsnwprintf_s_l(
                 wchar_t*       const _Buffer,
                                         size_t         const _BufferCount,
                                         size_t         const _MaxCount,
        wchar_t const* const _Format,
                                     _locale_t      const _Locale,
                                            va_list              _ArgList
    )



{
    int const _Result = __stdio_common_vsnwprintf_s(
        (*__local_stdio_printf_options()),
        _Buffer, _BufferCount, _MaxCount, _Format, _Locale, _ArgList);

    return _Result < 0 ? -1 : _Result;
}


 

__inline int __cdecl _vsnwprintf_s(
       wchar_t*       const _Buffer,
                               size_t         const _BufferCount,
                               size_t         const _MaxCount,
        wchar_t const* const _Format,
                                  va_list              _ArgList
    )



{
    return _vsnwprintf_s_l(_Buffer, _BufferCount, _MaxCount, _Format, 0, _ArgList);
}


__declspec(deprecated("This function or variable may be unsafe. Consider using " "_snwprintf_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __inline   int __cdecl _snwprintf(    wchar_t *_Buffer,   size_t _BufferCount,     wchar_t const* _Format, ...); __declspec(deprecated("This function or variable may be unsafe. Consider using " "_vsnwprintf_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __inline   int __cdecl _vsnwprintf(    wchar_t *_Buffer,   size_t _BufferCount,     wchar_t const* _Format, va_list _Args);

#pragma warning(pop)

 
 __declspec(deprecated("This function or variable may be unsafe. Consider using " "_vsnwprintf_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__inline int __cdecl _vsnwprintf(
        wchar_t*       _Buffer,
                                          size_t         _BufferCount,
                   wchar_t const* _Format,
                                             va_list        _ArgList
    )



{
    #pragma warning(push)
    #pragma warning(disable: 4996) 
    return _vsnwprintf_l(_Buffer, _BufferCount, _Format, 0, _ArgList);
    #pragma warning(pop)
}


extern "C++" { template <size_t _Size> inline   int __cdecl _vsnwprintf_s(  wchar_t (&_Buffer)[_Size],   size_t _BufferCount,     wchar_t const* _Format, va_list _ArgList) throw() { return _vsnwprintf_s(_Buffer, _Size, _BufferCount, _Format, _ArgList); } }

 

__inline int __cdecl _vswprintf_c_l(
                 wchar_t*       const _Buffer,
                                         size_t         const _BufferCount,
        wchar_t const* const _Format,
                                     _locale_t      const _Locale,
                                            va_list              _ArgList
    )



{
    int const _Result = __stdio_common_vswprintf(
        (*__local_stdio_printf_options()),
        _Buffer, _BufferCount, _Format, _Locale, _ArgList);

    return _Result < 0 ? -1 : _Result;
}


 

__inline int __cdecl _vswprintf_c(
       wchar_t*       const _Buffer,
                               size_t         const _BufferCount,
        wchar_t const* const _Format,
                                  va_list              _ArgList
    )



{
    return _vswprintf_c_l(_Buffer, _BufferCount, _Format, 0, _ArgList);
}


 

__inline int __cdecl _vswprintf_l(
                         wchar_t*       const _Buffer,
                                         size_t         const _BufferCount,
        wchar_t const* const _Format,
                                     _locale_t      const _Locale,
                                            va_list              _ArgList
    )



{
    #pragma warning(push)
    #pragma warning(disable: 4996) 
    return _vswprintf_c_l(_Buffer, _BufferCount, _Format, _Locale, _ArgList);
    #pragma warning(pop)
}


 

__inline int __cdecl __vswprintf_l(
                         wchar_t*       const _Buffer,
        wchar_t const* const _Format,
                                     _locale_t      const _Locale,
                                            va_list              _ArgList
    )



{
    return _vswprintf_l(_Buffer, (size_t)-1, _Format, _Locale, _ArgList);
}


 

__inline int __cdecl _vswprintf(
               wchar_t*       const _Buffer,
        wchar_t const* const _Format,
                                  va_list              _ArgList
    )



{
    return _vswprintf_l(_Buffer, (size_t)-1, _Format, 0, _ArgList);
}


 

__inline int __cdecl vswprintf(
                         wchar_t*       const _Buffer,
                                         size_t         const _BufferCount,
        wchar_t const* const _Format,
                                            va_list              _ArgList
    )



{
    return _vswprintf_c_l(_Buffer, _BufferCount, _Format, 0, _ArgList);
}


 

__inline int __cdecl _vswprintf_s_l(
                 wchar_t*       const _Buffer,
                                         size_t         const _BufferCount,
        wchar_t const* const _Format,
                                     _locale_t      const _Locale,
                                            va_list              _ArgList
    )



{
    int const _Result = __stdio_common_vswprintf_s(
        (*__local_stdio_printf_options()),
        _Buffer, _BufferCount, _Format, _Locale, _ArgList);

    return _Result < 0 ? -1 : _Result;
}




     
    __inline int __cdecl vswprintf_s(
           wchar_t*       const _Buffer,
                                   size_t         const _BufferCount,
            wchar_t const* const _Format,
                                      va_list              _ArgList
        )



    {
        return _vswprintf_s_l(_Buffer, _BufferCount, _Format, 0, _ArgList);
    }




extern "C++" { template <size_t _Size> inline   int __cdecl vswprintf_s(  wchar_t (&_Buffer)[_Size],     wchar_t const* _Format, va_list _ArgList) throw() { return vswprintf_s(_Buffer, _Size, _Format, _ArgList); } }

 

__inline int __cdecl _vswprintf_p_l(
                 wchar_t*       const _Buffer,
                                         size_t         const _BufferCount,
        wchar_t const* const _Format,
                                     _locale_t      const _Locale,
                                            va_list              _ArgList
    )



{
    int const _Result = __stdio_common_vswprintf_p(
        (*__local_stdio_printf_options()),
        _Buffer, _BufferCount, _Format, _Locale, _ArgList);

    return _Result < 0 ? -1 : _Result;
}


 

__inline int __cdecl _vswprintf_p(
       wchar_t*       const _Buffer,
                               size_t         const _BufferCount,
        wchar_t const* const _Format,
                                  va_list              _ArgList
    )



{
    return _vswprintf_p_l(_Buffer, _BufferCount, _Format, 0, _ArgList);
}


 
 
__inline int __cdecl _vscwprintf_l(
        wchar_t const* const _Format,
                                     _locale_t      const _Locale,
                                            va_list              _ArgList
    )



{
    int const _Result = __stdio_common_vswprintf(
        (*__local_stdio_printf_options()) | (1ULL << 1),
        0, 0, _Format, _Locale, _ArgList);

    return _Result < 0 ? -1 : _Result;
}


 
 
__inline int __cdecl _vscwprintf(
        wchar_t const* const _Format,
                                  va_list              _ArgList
    )



{
    return _vscwprintf_l(_Format, 0, _ArgList);
}


 
 
__inline int __cdecl _vscwprintf_p_l(
        wchar_t const* const _Format,
                                     _locale_t      const _Locale,
                                            va_list              _ArgList
    )



{
    int const _Result = __stdio_common_vswprintf_p(
        (*__local_stdio_printf_options()) | (1ULL << 1),
        0, 0, _Format, _Locale, _ArgList);

    return _Result < 0 ? -1 : _Result;
}


 
 
__inline int __cdecl _vscwprintf_p(
        wchar_t const* const _Format, 
                                  va_list              _ArgList
    )



{
    return _vscwprintf_p_l(_Format, 0, _ArgList);
}


 

__inline int __cdecl __swprintf_l(
                         wchar_t*       const _Buffer,
        wchar_t const* const _Format,
                                     _locale_t      const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = __vswprintf_l(_Buffer, _Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}


 

__inline int __cdecl _swprintf_l(
                         wchar_t*       const _Buffer,
                                         size_t         const _BufferCount,
        wchar_t const* const _Format,
                                     _locale_t      const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vswprintf_c_l(_Buffer, _BufferCount, _Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}


 

__inline int __cdecl _swprintf(
               wchar_t*       const _Buffer,
        wchar_t const* const _Format,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
    _Result = __vswprintf_l(_Buffer, _Format, 0, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}


 

__inline int __cdecl swprintf(
               wchar_t*       const _Buffer,
                               size_t         const _BufferCount,
        wchar_t const* const _Format,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
    _Result = _vswprintf_c_l(_Buffer, _BufferCount, _Format, 0, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}


#pragma warning(push)


#pragma warning(disable:4793 4996)

__declspec(deprecated("This function or variable may be unsafe. Consider using " "__swprintf_l_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __inline   int __cdecl __swprintf_l(    wchar_t *_Buffer,     wchar_t const* _Format,   _locale_t _Locale, ...); __declspec(deprecated("This function or variable may be unsafe. Consider using " "_vswprintf_s_l" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __inline   int __cdecl __vswprintf_l(    wchar_t *_Buffer,     wchar_t const* _Format,   _locale_t _Locale, va_list _Args);

__declspec(deprecated("This function or variable may be unsafe. Consider using " "_swprintf_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __inline   int __cdecl _swprintf(    wchar_t *_Buffer,     wchar_t const* _Format, ...); __declspec(deprecated("This function or variable may be unsafe. Consider using " "vswprintf_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __inline   int __cdecl _vswprintf(    wchar_t *_Buffer,     wchar_t const* _Format, va_list _Args);

#pragma warning(pop)

 

__inline int __cdecl _swprintf_s_l(
                 wchar_t*       const _Buffer,
                                         size_t         const _BufferCount,
        wchar_t const* const _Format,
                                     _locale_t      const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vswprintf_s_l(_Buffer, _BufferCount, _Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}




     
    __inline int __cdecl swprintf_s(
           wchar_t*       const _Buffer,
                                   size_t         const _BufferCount,
            wchar_t const* const _Format,
        ...)



    {
        int _Result;
        va_list _ArgList;
        ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
        _Result = _vswprintf_s_l(_Buffer, _BufferCount, _Format, 0, _ArgList);
        ((void)(_ArgList = (va_list)0));
        return _Result;
    }




extern "C++" { __pragma(warning(push)); __pragma(warning(disable: 4793)); template <size_t _Size> inline   int __cdecl swprintf_s(  wchar_t (&_Buffer)[_Size],     wchar_t const* _Format, ...) throw() { va_list _ArgList; ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format))))); return vswprintf_s(_Buffer, _Size, _Format, _ArgList); } __pragma(warning(pop)); }

 

__inline int __cdecl _swprintf_p_l(
                 wchar_t*       const _Buffer,
                                         size_t         const _BufferCount,
        wchar_t const* const _Format,
                                     _locale_t      const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vswprintf_p_l(_Buffer, _BufferCount, _Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}


 

__inline int __cdecl _swprintf_p(
       wchar_t*       const _Buffer,
                               size_t         const _BufferCount,
        wchar_t const* const _Format,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
    _Result = _vswprintf_p_l(_Buffer, _BufferCount, _Format, 0, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}


 

__inline int __cdecl _swprintf_c_l(
                 wchar_t*       const _Buffer,
                                         size_t         const _BufferCount,
        wchar_t const* const _Format,
                                     _locale_t      const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vswprintf_c_l(_Buffer, _BufferCount, _Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}


 

__inline int __cdecl _swprintf_c(
       wchar_t*       const _Buffer,
                               size_t         const _BufferCount,
        wchar_t const* const _Format,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
    _Result = _vswprintf_c_l(_Buffer, _BufferCount, _Format, 0, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}


 
 __declspec(deprecated("This function or variable may be unsafe. Consider using " "_snwprintf_s_l" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__inline int __cdecl _snwprintf_l(
        wchar_t*       const _Buffer,
                                          size_t         const _BufferCount,
         wchar_t const* const _Format,
                                      _locale_t      const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));

    #pragma warning(push)
    #pragma warning(disable: 4996) 
    _Result = _vsnwprintf_l(_Buffer, _BufferCount, _Format, _Locale, _ArgList);
    #pragma warning(pop)

    ((void)(_ArgList = (va_list)0));
    return _Result;
}


 

__inline int __cdecl _snwprintf(
        wchar_t*       _Buffer,
                                          size_t         _BufferCount,
                   wchar_t const* _Format,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));

    #pragma warning(push)
    #pragma warning(disable: 4996) 
    _Result = _vsnwprintf_l(_Buffer, _BufferCount, _Format, 0, _ArgList);
    #pragma warning(pop)

    ((void)(_ArgList = (va_list)0));
    return _Result;
}


 

__inline int __cdecl _snwprintf_s_l(
                 wchar_t*       const _Buffer,
                                         size_t         const _BufferCount,
                                         size_t         const _MaxCount,
        wchar_t const* const _Format,
                                     _locale_t      const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vsnwprintf_s_l(_Buffer, _BufferCount, _MaxCount, _Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}


 

__inline int __cdecl _snwprintf_s(
       wchar_t*       const _Buffer,
                               size_t         const _BufferCount,
                               size_t         const _MaxCount,
        wchar_t const* const _Format,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
    _Result = _vsnwprintf_s_l(_Buffer, _BufferCount, _MaxCount, _Format, 0, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}


extern "C++" { __pragma(warning(push)); __pragma(warning(disable: 4793)); template <size_t _Size> inline   int __cdecl _snwprintf_s(  wchar_t (&_Buffer)[_Size],   size_t _BufferCount,     wchar_t const* _Format, ...) throw() { va_list _ArgList; ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format))))); return _vsnwprintf_s(_Buffer, _Size, _BufferCount, _Format, _ArgList); } __pragma(warning(pop)); }

 
__inline int __cdecl _scwprintf_l(
        wchar_t const* const _Format,
                                     _locale_t      const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vscwprintf_l(_Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}


 
 
__inline int __cdecl _scwprintf(
        wchar_t const* const _Format,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
    _Result = _vscwprintf_l(_Format, 0, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}


 
 
__inline int __cdecl _scwprintf_p_l(
        wchar_t const* const _Format,
                                     _locale_t      const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vscwprintf_p_l(_Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}


 
 
__inline int __cdecl _scwprintf_p(
        wchar_t const* const _Format,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
    _Result = _vscwprintf_p_l(_Format, 0, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}




    #pragma warning(push)
    #pragma warning(disable: 4141 4412 4793 4996 6054)

    

        extern "C++" __declspec(deprecated("function has been changed to conform with the ISO C standard, " "adding an extra character count parameter. To use the traditional " "Microsoft version, set _CRT_NON_CONFORMING_SWPRINTFS.")) __declspec(deprecated("This function or variable may be unsafe. Consider using " "swprintf_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
        inline int swprintf(
                       wchar_t*       const _Buffer,
                wchar_t const* const _Format,
            ...) throw()
        {
            int _Result;
            va_list _ArgList;
            ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
            #pragma warning(suppress: 28719)
            _Result = vswprintf(_Buffer, 2147483647, _Format, _ArgList);       
            ((void)(_ArgList = (va_list)0));
            return _Result;
        }

        extern "C++" __declspec(deprecated("function has been changed to conform with the ISO C standard, " "adding an extra character count parameter. To use the traditional " "Microsoft version, set _CRT_NON_CONFORMING_SWPRINTFS.")) __declspec(deprecated("This function or variable may be unsafe. Consider using " "vswprintf_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
        inline int __cdecl vswprintf(
                       wchar_t*       const _Buffer,
                wchar_t const* const _Format,
                                          va_list              _ArgList
            ) throw()
        {
            #pragma warning(suppress: 28719)
            return vswprintf(_Buffer, 2147483647, _Format, _ArgList);
        }

        extern "C++" __declspec(deprecated("function has been changed to conform with the ISO C standard, " "adding an extra character count parameter. To use the traditional " "Microsoft version, set _CRT_NON_CONFORMING_SWPRINTFS.")) __declspec(deprecated("This function or variable may be unsafe. Consider using " "_swprintf_s_l" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
        inline int _swprintf_l(
                                 wchar_t*       const _Buffer,
                wchar_t const* const _Format,
                                             _locale_t      const _Locale,
            ...) throw()
        {
            int _Result;
            va_list _ArgList;
            ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
            _Result = _vswprintf_l(_Buffer, (size_t)-1, _Format, _Locale, _ArgList);
            ((void)(_ArgList = (va_list)0));
            return _Result;
        }

        extern "C++" __declspec(deprecated("function has been changed to conform with the ISO C standard, " "adding an extra character count parameter. To use the traditional " "Microsoft version, set _CRT_NON_CONFORMING_SWPRINTFS.")) __declspec(deprecated("This function or variable may be unsafe. Consider using " "_vswprintf_s_l" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
        inline int __cdecl _vswprintf_l(
                                 wchar_t*       const _Buffer,
                wchar_t const* const _Format,
                                             _locale_t      const _Locale,
                                                    va_list              _ArgList
            ) throw()
        {
            return _vswprintf_l(_Buffer, (size_t)-1, _Format, _Locale, _ArgList);
        }

    

    #pragma warning(pop)















 
__declspec(dllimport) int __cdecl __stdio_common_vswscanf(
                                        unsigned __int64 _Options,
              wchar_t const*   _Buffer,
                                        size_t           _BufferCount,
        wchar_t const*   _Format,
                                    _locale_t        _Locale,
                                           va_list          _ArgList
    );

 

__inline int __cdecl _vswscanf_l(
                             wchar_t const* const _Buffer,
        wchar_t const* const _Format,
                           _locale_t      const _Locale,
                                  va_list              _ArgList
    )



{
    return __stdio_common_vswscanf(
        (*__local_stdio_scanf_options ()),
        _Buffer, (size_t)-1, _Format, _Locale, _ArgList);
}


 

__inline int __cdecl vswscanf(
                             wchar_t const* _Buffer,
        wchar_t const* _Format,
                                  va_list        _ArgList
    )



{
    return _vswscanf_l(_Buffer, _Format, 0, _ArgList);
}


 

__inline int __cdecl _vswscanf_s_l(
                             wchar_t const* const _Buffer,
        wchar_t const* const _Format,
                           _locale_t      const _Locale,
                                  va_list              _ArgList
    )



{
    return __stdio_common_vswscanf(
        (*__local_stdio_scanf_options ()) | (1ULL << 0),
        _Buffer, (size_t)-1, _Format, _Locale, _ArgList);
}




     
    
    __inline int __cdecl vswscanf_s(
                                 wchar_t const* const _Buffer,
            wchar_t const* const _Format,
                                      va_list              _ArgList
        )



    {
        return _vswscanf_s_l(_Buffer, _Format, 0, _ArgList);
    }




extern "C++" { template <size_t _Size> inline   int __cdecl vswscanf_s(  wchar_t (&_Buffer)[_Size],     wchar_t const* _Format, va_list _Args) throw() { return vswscanf_s(_Buffer, _Size, _Format, _Args); } }

 
 __declspec(deprecated("This function or variable may be unsafe. Consider using " "_vsnwscanf_s_l" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__inline int __cdecl _vsnwscanf_l(
              wchar_t const* const _Buffer,
                                        size_t         const _BufferCount,
        wchar_t const* const _Format,
                                    _locale_t      const _Locale,
                                           va_list              _ArgList
    )



{
    return __stdio_common_vswscanf(
        (*__local_stdio_scanf_options ()),
        _Buffer, _BufferCount, _Format, _Locale, _ArgList);
}


 

__inline int __cdecl _vsnwscanf_s_l(
                wchar_t const* const _Buffer,
                                          size_t         const _BufferCount,
        wchar_t const* const _Format,
                                      _locale_t      const _Locale,
                                             va_list              _ArgList
    )



{
    return __stdio_common_vswscanf(
        (*__local_stdio_scanf_options ()) | (1ULL << 0),
        _Buffer, _BufferCount, _Format, _Locale, _ArgList);
}


 
 __declspec(deprecated("This function or variable may be unsafe. Consider using " "_swscanf_s_l" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__inline int __cdecl _swscanf_l(
                                      wchar_t const* const _Buffer,
        wchar_t const* const _Format,
                                    _locale_t            _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vswscanf_l(_Buffer, _Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}


 
  __declspec(deprecated("This function or variable may be unsafe. Consider using " "swscanf_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__inline int __cdecl swscanf(
                            wchar_t const* const _Buffer,
        wchar_t const* const _Format,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
    _Result = _vswscanf_l(_Buffer, _Format, 0, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}


 

__inline int __cdecl _swscanf_s_l(
                                        wchar_t const* const _Buffer,
        wchar_t const* const _Format,
                                      _locale_t      const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vswscanf_s_l(_Buffer, _Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}




     
    
    __inline int __cdecl swscanf_s(
                                  wchar_t const* const _Buffer,
            wchar_t const* const _Format,
        ...)



    {
        int _Result;
        va_list _ArgList;
        ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
        _Result = _vswscanf_s_l(_Buffer, _Format, 0, _ArgList);  
        ((void)(_ArgList = (va_list)0));
        return _Result;
    }




 
 __declspec(deprecated("This function or variable may be unsafe. Consider using " "_snwscanf_s_l" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__inline int __cdecl _snwscanf_l(
              wchar_t const* const _Buffer,
                                        size_t         const _BufferCount,
        wchar_t const* const _Format,
                                    _locale_t      const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));

    #pragma warning(push)
    #pragma warning(disable: 4996) 
    _Result = _vsnwscanf_l(_Buffer, _BufferCount, _Format, _Locale, _ArgList);
    #pragma warning(pop)

    ((void)(_ArgList = (va_list)0));
    return _Result;
}


 
 __declspec(deprecated("This function or variable may be unsafe. Consider using " "_snwscanf_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__inline int __cdecl _snwscanf(
        wchar_t const* const _Buffer,
                                  size_t         const _BufferCount,
            wchar_t const* const _Format,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));

    #pragma warning(push)
    #pragma warning(disable: 4996) 
    _Result = _vsnwscanf_l(_Buffer, _BufferCount, _Format, 0, _ArgList);
    #pragma warning(pop)

    ((void)(_ArgList = (va_list)0));
    return _Result;
}


 

__inline int __cdecl _snwscanf_s_l(
                wchar_t const* const _Buffer,
                                          size_t         const _BufferCount,
        wchar_t const* const _Format,
                                      _locale_t      const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vsnwscanf_s_l(_Buffer, _BufferCount, _Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}


 

__inline int __cdecl _snwscanf_s(
         wchar_t const* const _Buffer,
                                   size_t         const _BufferCount,
           wchar_t const* const _Format,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
    _Result = _vsnwscanf_s_l(_Buffer, _BufferCount, _Format, 0, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}








} __pragma(pack(pop))












#pragma once






__pragma(pack(push, 8)) extern "C" {










    
    __declspec(dllimport) errno_t __cdecl wcscat_s(
          wchar_t* _Destination,
          rsize_t _SizeInWords,
          wchar_t const* _Source
        );

    
    __declspec(dllimport) errno_t __cdecl wcscpy_s(
          wchar_t* _Destination,
          rsize_t _SizeInWords,
          wchar_t const* _Source
        );
    
    
    __declspec(dllimport) errno_t __cdecl wcsncat_s(
          wchar_t*       _Destination,
                                     rsize_t        _SizeInWords,
               wchar_t const* _Source,
                                     rsize_t        _MaxCount
        );
    
    
    __declspec(dllimport) errno_t __cdecl wcsncpy_s(
          wchar_t*       _Destination,
                                  rsize_t        _SizeInWords,
            wchar_t const* _Source,
                                  rsize_t        _MaxCount
        );
    
     
    __declspec(dllimport) wchar_t* __cdecl wcstok_s(
                          wchar_t*       _String,
                                 wchar_t const* _Delimiter,
            wchar_t**      _Context
        );















 
__declspec(dllimport) __declspec(allocator) wchar_t* __cdecl _wcsdup(
      wchar_t const* _String
    );







extern "C++" { template <size_t _Size> inline errno_t __cdecl wcscat_s(wchar_t (&_Destination)[_Size],   wchar_t const* _Source) throw() { return wcscat_s(_Destination, _Size, _Source); } }



    __declspec(deprecated("This function or variable may be unsafe. Consider using " "wcscat_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport) wchar_t* __cdecl wcscat( wchar_t *_Destination,  wchar_t const* _Source);



 
__declspec(dllimport) int __cdecl wcscmp(
      wchar_t const* _String1,
      wchar_t const* _String2
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl wcscpy_s(wchar_t (&_Destination)[_Size],   wchar_t const* _Source) throw() { return wcscpy_s(_Destination, _Size, _Source); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "wcscpy_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport) wchar_t* __cdecl wcscpy( wchar_t *_Destination,  wchar_t const* _Source);

 
__declspec(dllimport) size_t __cdecl wcscspn(
      wchar_t const* _String,
      wchar_t const* _Control
    );

 
__declspec(dllimport) size_t __cdecl wcslen(
      wchar_t const* _String
    );

 


__declspec(dllimport) size_t __cdecl wcsnlen(
      wchar_t const* _Source,
                            size_t         _MaxCount
    );



     
    
    
    static __inline size_t __cdecl wcsnlen_s(
          wchar_t const* _Source,
                                size_t         _MaxCount
        )
    {
        return (_Source == 0) ? 0 : wcsnlen(_Source, _MaxCount);
    }



extern "C++" { template <size_t _Size> inline errno_t __cdecl wcsncat_s(  wchar_t (&_Destination)[_Size],   wchar_t const* _Source,   size_t _Count) throw() { return wcsncat_s(_Destination, _Size, _Source, _Count); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "wcsncat_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport) wchar_t* __cdecl wcsncat(  wchar_t *_Destination,   wchar_t const* _Source,   size_t _Count);

 
__declspec(dllimport) int __cdecl wcsncmp(
      wchar_t const* _String1,
      wchar_t const* _String2,
                            size_t         _MaxCount
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl wcsncpy_s(wchar_t (&_Destination)[_Size],   wchar_t const* _Source,   size_t _Count) throw() { return wcsncpy_s(_Destination, _Size, _Source, _Count); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "wcsncpy_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport) wchar_t* __cdecl wcsncpy(    wchar_t *_Destination,   wchar_t const* _Source,   size_t _Count);

 
__declspec(dllimport) wchar_t const* __cdecl wcspbrk(
      wchar_t const* _String,
      wchar_t const* _Control
    );

 
__declspec(dllimport) size_t __cdecl wcsspn(
      wchar_t const* _String,
      wchar_t const* _Control
    );

  __declspec(deprecated("This function or variable may be unsafe. Consider using " "wcstok_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__declspec(dllimport) wchar_t* __cdecl wcstok(
                          wchar_t*       _String,
                                 wchar_t const* _Delimiter,
        wchar_t**      _Context
    );



    


        



    

    #pragma warning(push)
    #pragma warning(disable: 4141 4996) 

      __declspec(deprecated("This function or variable may be unsafe. Consider using " "wcstok_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
    static __inline wchar_t* __cdecl _wcstok(
          wchar_t*       const _String,
                 wchar_t const* const _Delimiter
        )
    {
        return wcstok(_String, _Delimiter, 0);
    }

    



    
        extern "C++"   __declspec(deprecated("wcstok has been changed to conform with the ISO C standard, " "adding an extra context parameter. To use the legacy Microsoft " "wcstok, define _CRT_NON_CONFORMING_WCSTOK.")) 
        inline wchar_t* __cdecl wcstok(
              wchar_t*       _String,
                     wchar_t const* _Delimiter
            ) throw()
        {
            return wcstok(_String, _Delimiter, 0);
        }
    

    #pragma warning(pop)





 
  __declspec(deprecated("This function or variable may be unsafe. Consider using " "_wcserror_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__declspec(dllimport) wchar_t* __cdecl _wcserror(
      int _ErrorNumber
    );


__declspec(dllimport) errno_t __cdecl _wcserror_s(
      wchar_t* _Buffer,
                                  size_t   _SizeInWords,
                                  int      _ErrorNumber
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl _wcserror_s(wchar_t (&_Buffer)[_Size],   int _Error) throw() { return _wcserror_s(_Buffer, _Size, _Error); } }

 
 
  __declspec(deprecated("This function or variable may be unsafe. Consider using " "__wcserror_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__declspec(dllimport) wchar_t* __cdecl __wcserror(
      wchar_t const* _String
    );

 __declspec(dllimport) errno_t __cdecl __wcserror_s(
      wchar_t*       _Buffer,
                                  size_t         _SizeInWords,
                                wchar_t const* _ErrorMessage
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl __wcserror_s(wchar_t (&_Buffer)[_Size],   wchar_t const* _ErrorMessage) throw() { return __wcserror_s(_Buffer, _Size, _ErrorMessage); } }

  __declspec(dllimport) int __cdecl _wcsicmp(
      wchar_t const* _String1,
      wchar_t const* _String2
    );

  __declspec(dllimport) int __cdecl _wcsicmp_l(
        wchar_t const* _String1,
        wchar_t const* _String2,
      _locale_t      _Locale
    );

  __declspec(dllimport) int __cdecl _wcsnicmp(
      wchar_t const* _String1,
      wchar_t const* _String2,
                            size_t         _MaxCount
    );

  __declspec(dllimport) int __cdecl _wcsnicmp_l(
      wchar_t const* _String1,
      wchar_t const* _String2,
                            size_t         _MaxCount,
                        _locale_t      _Locale
    );

 __declspec(dllimport) errno_t __cdecl _wcsnset_s(
      wchar_t* _Destination,
                                 size_t   _SizeInWords,
                                 wchar_t  _Value,
                                 size_t   _MaxCount
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl _wcsnset_s(  wchar_t (&_Destination)[_Size],   wchar_t _Value,   size_t _MaxCount) throw() { return _wcsnset_s(_Destination, _Size, _Value, _MaxCount); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "_wcsnset_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport) wchar_t* __cdecl _wcsnset(  wchar_t *_String,   wchar_t _Value,   size_t _MaxCount);

__declspec(dllimport) wchar_t* __cdecl _wcsrev(
      wchar_t* _String
    );

 __declspec(dllimport) errno_t __cdecl _wcsset_s(
      wchar_t* _Destination,
                                 size_t   _SizeInWords,
                                 wchar_t  _Value
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl _wcsset_s(  wchar_t (&_String)[_Size],   wchar_t _Value) throw() { return _wcsset_s(_String, _Size, _Value); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "_wcsset_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport) wchar_t* __cdecl _wcsset(  wchar_t *_String,   wchar_t _Value);

 __declspec(dllimport) errno_t __cdecl _wcslwr_s(
      wchar_t* _String,
                                 size_t   _SizeInWords
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl _wcslwr_s(  wchar_t (&_String)[_Size]) throw() { return _wcslwr_s(_String, _Size); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "_wcslwr_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport) wchar_t* __cdecl _wcslwr( wchar_t *_String);


__declspec(dllimport) errno_t __cdecl _wcslwr_s_l(
      wchar_t*  _String,
                                 size_t    _SizeInWords,
                             _locale_t _Locale
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl _wcslwr_s_l(  wchar_t (&_String)[_Size],   _locale_t _Locale) throw() { return _wcslwr_s_l(_String, _Size, _Locale); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "_wcslwr_s_l" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport) wchar_t* __cdecl _wcslwr_l(  wchar_t *_String,   _locale_t _Locale);


__declspec(dllimport) errno_t __cdecl _wcsupr_s(
      wchar_t* _String,
                          size_t   _Size
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl _wcsupr_s(  wchar_t (&_String)[_Size]) throw() { return _wcsupr_s(_String, _Size); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "_wcsupr_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport) wchar_t* __cdecl _wcsupr( wchar_t *_String);


__declspec(dllimport) errno_t __cdecl _wcsupr_s_l(
      wchar_t*  _String,
                          size_t    _Size,
                      _locale_t _Locale
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl _wcsupr_s_l(  wchar_t (&_String)[_Size],   _locale_t _Locale) throw() { return _wcsupr_s_l(_String, _Size, _Locale); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "_wcsupr_s_l" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport) wchar_t* __cdecl _wcsupr_l(  wchar_t *_String,   _locale_t _Locale);

 

__declspec(dllimport) size_t __cdecl wcsxfrm(
        wchar_t*       _Destination,
                                         wchar_t const* _Source,
                size_t         _MaxCount
    );

 

__declspec(dllimport) size_t __cdecl _wcsxfrm_l(
        wchar_t*       _Destination,
                                         wchar_t const* _Source,
                size_t         _MaxCount,
                                       _locale_t      _Locale
    );

 
__declspec(dllimport) int __cdecl wcscoll(
      wchar_t const* _String1,
      wchar_t const* _String2
    );

 
__declspec(dllimport) int __cdecl _wcscoll_l(
        wchar_t const* _String1,
        wchar_t const* _String2,
      _locale_t      _Locale
    );

 
__declspec(dllimport) int __cdecl _wcsicoll(
      wchar_t const* _String1,
      wchar_t const* _String2
    );

 
__declspec(dllimport) int __cdecl _wcsicoll_l(
        wchar_t const* _String1,
        wchar_t const* _String2,
      _locale_t      _Locale
    );

 
__declspec(dllimport) int __cdecl _wcsncoll(
      wchar_t const* _String1,
      wchar_t const* _String2,
                            size_t         _MaxCount
    );

 
__declspec(dllimport) int __cdecl _wcsncoll_l(
      wchar_t const* _String1,
      wchar_t const* _String2,
                            size_t         _MaxCount,
                        _locale_t      _Locale
    );

 
__declspec(dllimport) int __cdecl _wcsnicoll(
      wchar_t const* _String1,
      wchar_t const* _String2,
                            size_t         _MaxCount
    );

 
__declspec(dllimport) int __cdecl _wcsnicoll_l(
      wchar_t const* _String1,
      wchar_t const* _String2,
                            size_t         _MaxCount,
                        _locale_t      _Locale
    );









extern "C++" {

     
    
    inline wchar_t* __cdecl wcschr(  wchar_t* _String, wchar_t _C)
    {
        return const_cast<wchar_t*>(wcschr(static_cast<wchar_t const*>(_String), _C));
    }

     
    inline wchar_t* __cdecl wcspbrk(  wchar_t* _String,   wchar_t const* _Control)
    {
        return const_cast<wchar_t*>(wcspbrk(static_cast<wchar_t const*>(_String), _Control));
    }

     
    inline wchar_t* __cdecl wcsrchr(  wchar_t* _String,   wchar_t _C)
    {
        return const_cast<wchar_t*>(wcsrchr(static_cast<wchar_t const*>(_String), _C));
    }

       
    
    inline wchar_t* __cdecl wcsstr(  wchar_t* _String,   wchar_t const*_SubStr)
    {
        return const_cast<wchar_t*>(wcsstr(static_cast<wchar_t const*>(_String), _SubStr));
    }

}










    




      __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_wcsdup" ". See online help for details."))
    __declspec(dllimport) wchar_t* __cdecl wcsdup(
          wchar_t const* _String
        );

    



    
    

      __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_wcsicmp" ". See online help for details."))
    __declspec(dllimport) int __cdecl wcsicmp(
          wchar_t const* _String1,
          wchar_t const* _String2
        );
    
      __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_wcsnicmp" ". See online help for details."))
    __declspec(dllimport) int __cdecl wcsnicmp(
          wchar_t const* _String1,
          wchar_t const* _String2,
                                size_t         _MaxCount
        );
    
    __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_wcsnset" ". See online help for details."))
     
    __declspec(dllimport) wchar_t* __cdecl wcsnset(
          wchar_t* _String,
                                  wchar_t  _Value,
                                  size_t   _MaxCount
        );
    
    __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_wcsrev" ". See online help for details."))
     
    __declspec(dllimport) wchar_t* __cdecl wcsrev(
          wchar_t* _String
        );
    
    __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_wcsset" ". See online help for details."))
     
    __declspec(dllimport) wchar_t* __cdecl wcsset(
          wchar_t* _String,
               wchar_t  _Value
        );
    
    __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_wcslwr" ". See online help for details."))
     
    __declspec(dllimport) wchar_t* __cdecl wcslwr(
          wchar_t* _String
        );
    
    __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_wcsupr" ". See online help for details."))
     
    __declspec(dllimport) wchar_t* __cdecl wcsupr(
          wchar_t* _String
        );
    
      __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_wcsicoll" ". See online help for details."))
    __declspec(dllimport) int __cdecl wcsicoll(
          wchar_t const* _String1,
          wchar_t const* _String2
        );





} __pragma(pack(pop))












#pragma once



__pragma(pack(push, 8)) extern "C" {








struct tm
{
    int tm_sec;   
    int tm_min;   
    int tm_hour;  
    int tm_mday;  
    int tm_mon;   
    int tm_year;  
    int tm_wday;  
    int tm_yday;  
    int tm_isdst; 
};







  __declspec(deprecated("This function or variable may be unsafe. Consider using " "_wasctime_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
 
 
__declspec(dllimport) wchar_t* __cdecl _wasctime(
      struct tm const* _Tm
    );

 

__declspec(dllimport) errno_t __cdecl _wasctime_s(
        wchar_t*         _Buffer,
                                          size_t           _SizeInWords,
                                                       struct tm const* _Tm
    );

extern "C++" { template <size_t _Size> inline   errno_t __cdecl _wasctime_s(  wchar_t (&_Buffer)[_Size],   struct tm const* _Time) throw() { return _wasctime_s(_Buffer, _Size, _Time); } }

 

__declspec(dllimport) size_t __cdecl wcsftime(
       wchar_t*         _Buffer,
                               size_t           _SizeInWords,
                             wchar_t const*   _Format,
                               struct tm const* _Tm
    );

 

__declspec(dllimport) size_t __cdecl _wcsftime_l(
       wchar_t*         _Buffer,
                               size_t           _SizeInWords,
                             wchar_t const*   _Format,
                               struct tm const* _Tm,
                           _locale_t        _Locale
    );

 
  __declspec(deprecated("This function or variable may be unsafe. Consider using " "_wctime32_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__declspec(dllimport) wchar_t* __cdecl _wctime32(
      __time32_t const* _Time
    );


__declspec(dllimport) errno_t __cdecl _wctime32_s(
        wchar_t*          _Buffer,
                                      size_t            _SizeInWords,
                                                       __time32_t const* _Time
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl _wctime32_s(  wchar_t (&_Buffer)[_Size],   __time32_t const* _Time) throw() { return _wctime32_s(_Buffer, _Size, _Time); } }

 
 
  __declspec(deprecated("This function or variable may be unsafe. Consider using " "_wctime64_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__declspec(dllimport) wchar_t* __cdecl _wctime64(
      __time64_t const* _Time
    );


__declspec(dllimport) errno_t __cdecl _wctime64_s(
        wchar_t*          _Buffer,
                                      size_t            _SizeInWords,
                                                       __time64_t const* _Time);

extern "C++" { template <size_t _Size> inline errno_t __cdecl _wctime64_s(  wchar_t (&_Buffer)[_Size],   __time64_t const* _Time) throw() { return _wctime64_s(_Buffer, _Size, _Time); } }


__declspec(dllimport) errno_t __cdecl _wstrdate_s(
        wchar_t* _Buffer,
                                      size_t   _SizeInWords
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl _wstrdate_s(  wchar_t (&_Buffer)[_Size]) throw() { return _wstrdate_s(_Buffer, _Size); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "_wstrdate_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport)  wchar_t* __cdecl _wstrdate( wchar_t *_Buffer);


__declspec(dllimport) errno_t __cdecl _wstrtime_s(
        wchar_t* _Buffer,
                                      size_t   _SizeInWords
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl _wstrtime_s(  wchar_t (&_Buffer)[_Size]) throw() { return _wstrtime_s(_Buffer, _Size); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "_wstrtime_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport)  wchar_t* __cdecl _wstrtime( wchar_t *_Buffer);









    #pragma warning(push)
    #pragma warning(disable: 4996)

    




















         
        static __inline wchar_t * __cdecl _wctime(
              time_t const* const _Time)
        {
            return _wctime64(_Time);
        }

        
        static __inline errno_t __cdecl _wctime_s(
                  wchar_t*      const _Buffer,
                                                             size_t        const _SizeInWords,
                                                             time_t const* const _Time
            )
        {
            return _wctime64_s(_Buffer, _SizeInWords, _Time);
        }

    

    #pragma warning(pop)




} __pragma(pack(pop))









#pragma once










#pragma once




    

    typedef unsigned short _ino_t; 

    
        typedef _ino_t ino_t;
    





    

    typedef unsigned int _dev_t; 

    
        typedef _dev_t dev_t;
    





    

    typedef long _off_t; 

    
        typedef _off_t off_t;
    



__pragma(pack(push, 8)) extern "C" {



#pragma warning(push)
#pragma warning(disable:4820) 








struct _stat32
{
    _dev_t         st_dev;
    _ino_t         st_ino;
    unsigned short st_mode;
    short          st_nlink;
    short          st_uid;
    short          st_gid;
    _dev_t         st_rdev;
    _off_t         st_size;
    __time32_t     st_atime;
    __time32_t     st_mtime;
    __time32_t     st_ctime;
};

struct _stat32i64
{
    _dev_t         st_dev;
    _ino_t         st_ino;
    unsigned short st_mode;
    short          st_nlink;
    short          st_uid;
    short          st_gid;
    _dev_t         st_rdev;
    __int64        st_size;
    __time32_t     st_atime;
    __time32_t     st_mtime;
    __time32_t     st_ctime;
};

struct _stat64i32
{
    _dev_t         st_dev;
    _ino_t         st_ino;
    unsigned short st_mode;
    short          st_nlink;
    short          st_uid;
    short          st_gid;
    _dev_t         st_rdev;
    _off_t         st_size;
    __time64_t     st_atime;
    __time64_t     st_mtime;
    __time64_t     st_ctime;
};

struct _stat64
{
    _dev_t         st_dev;
    _ino_t         st_ino;
    unsigned short st_mode;
    short          st_nlink;
    short          st_uid;
    short          st_gid;
    _dev_t         st_rdev;
    __int64        st_size;
    __time64_t     st_atime;
    __time64_t     st_mtime;
    __time64_t     st_ctime;
};




    struct stat
    {
        _dev_t         st_dev;
        _ino_t         st_ino;
        unsigned short st_mode;
        short          st_nlink;
        short          st_uid;
        short          st_gid;
        _dev_t         st_rdev;
        _off_t         st_size;
        time_t         st_atime;
        time_t         st_mtime;
        time_t         st_ctime;
    };



















    
    
    
    
    
    
    

















    
    
    
    
    
    




__declspec(dllimport) int __cdecl _fstat32(
       int             _FileHandle,
      struct _stat32* _Stat
    );

__declspec(dllimport) int __cdecl _fstat32i64(
       int                _FileHandle,
      struct _stat32i64* _Stat
    );

__declspec(dllimport) int __cdecl _fstat64i32(
       int                _FileHandle,
      struct _stat64i32* _Stat
    );

__declspec(dllimport) int __cdecl _fstat64(
       int             _FileHandle,
      struct _stat64* _Stat
    );

__declspec(dllimport) int __cdecl _stat32(
      char const*     _FileName,
       struct _stat32* _Stat
    );

__declspec(dllimport) int __cdecl _stat32i64(
      char const*        _FileName,
       struct _stat32i64* _Stat
    );

__declspec(dllimport) int __cdecl _stat64i32(
      char const*        _FileName,
       struct _stat64i32* _Stat
    );

__declspec(dllimport) int __cdecl _stat64(
      char const*     _FileName,
       struct _stat64* _Stat
    );

__declspec(dllimport) int __cdecl _wstat32(
      wchar_t const*  _FileName,
       struct _stat32* _Stat
    );

__declspec(dllimport) int __cdecl _wstat32i64(
      wchar_t const*     _FileName,
       struct _stat32i64* _Stat
    );

__declspec(dllimport) int __cdecl _wstat64i32(
      wchar_t const*     _FileName,
       struct _stat64i32* _Stat
    );

__declspec(dllimport) int __cdecl _wstat64(
      wchar_t const*  _FileName,
       struct _stat64* _Stat
    );




    















        static __inline int __cdecl fstat(int const _FileHandle, struct stat* const _Stat)
        {
            typedef char __static_assert_t[(sizeof(struct stat) == sizeof(struct _stat64i32)) != 0];
            return _fstat64i32(_FileHandle, (struct _stat64i32*)_Stat);
        }
        static __inline int __cdecl stat(char const* const _FileName, struct stat* const _Stat)
        {
            typedef char __static_assert_t[(sizeof(struct stat) == sizeof(struct _stat64i32)) != 0];
            return _stat64i32(_FileName, (struct _stat64i32*)_Stat);
        }

    




#pragma warning(pop)



} __pragma(pack(pop))




__pragma(pack(push, 8)) extern "C" {








typedef wchar_t _Wint_t;




__declspec(dllimport) wchar_t* __cdecl _wsetlocale(
            int            _Category,
      wchar_t const* _Locale
    );


__declspec(dllimport) _locale_t __cdecl _wcreate_locale(
        int            _Category,
      wchar_t const* _Locale
    );



__declspec(dllimport) wint_t __cdecl btowc(
      int _Ch
    );

__declspec(dllimport) size_t __cdecl mbrlen(
        char const* _Ch,
                                                size_t      _SizeInBytes,
                                             mbstate_t*  _State
    );

__declspec(dllimport) size_t __cdecl mbrtowc(
                              wchar_t*    _DstCh,
        char const* _SrcCh,
                                                size_t      _SizeInBytes,
                                             mbstate_t*  _State
    );

 
__declspec(dllimport) errno_t __cdecl mbsrtowcs_s(
                              size_t*      _Retval,
              wchar_t*     _Dst,
                                   size_t       _Size,
                      char const** _PSrc,
                                   size_t       _N,
                                mbstate_t*   _State
    );

extern "C++" { template <size_t _Size> inline   errno_t __cdecl mbsrtowcs_s(  size_t* _Retval,   wchar_t (&_Dest)[_Size],     char const** _PSource,   size_t _Count,   mbstate_t* _State) throw() { return mbsrtowcs_s(_Retval, _Dest, _Size, _PSource, _Count, _State); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "mbsrtowcs_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))  __declspec(dllimport) size_t __cdecl mbsrtowcs( wchar_t *_Dest,  char const** _PSrc,  size_t _Count,  mbstate_t* _State);

 
__declspec(dllimport) errno_t __cdecl wcrtomb_s(
                             size_t*    _Retval,
      char*      _Dst,
                                  size_t     _SizeInBytes,
                                  wchar_t    _Ch,
                           mbstate_t* _State
    );

extern "C++" { template <size_t _Size> inline   errno_t __cdecl wcrtomb_s(  size_t* _Retval,   char (&_Dest)[_Size],   wchar_t _Source,   mbstate_t* _State) throw() { return wcrtomb_s(_Retval, _Dest, _Size, _Source, _State); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "wcrtomb_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport) size_t __cdecl wcrtomb(  char *_Dest,  wchar_t _Source,  mbstate_t* _State);

 
__declspec(dllimport) errno_t __cdecl wcsrtombs_s(
                                              size_t*         _Retval,
      char*           _Dst,
                                                   size_t          _SizeInBytes,
                                wchar_t const** _Src,
                                                   size_t          _Size,
                                            mbstate_t*      _State
    );

extern "C++" { template <size_t _Size> inline   errno_t __cdecl wcsrtombs_s(  size_t* _Retval,   char (&_Dest)[_Size],     wchar_t const** _PSrc,   size_t _Count,   mbstate_t* _State) throw() { return wcsrtombs_s(_Retval, _Dest, _Size, _PSrc, _Count, _State); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "wcsrtombs_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport) size_t __cdecl wcsrtombs(  char *_Dest,  wchar_t const** _PSource,  size_t _Count,  mbstate_t* _State);

__declspec(dllimport) int __cdecl wctob(
      wint_t _WCh
    );





    

         
        errno_t __cdecl wmemcpy_s(
              wchar_t*       _S1,
                                      rsize_t        _N1,
                        wchar_t const* _S2,
                                      rsize_t        _N
            );

         
        errno_t __cdecl wmemmove_s(
              wchar_t*       _S1,
                                      rsize_t        _N1,
                        wchar_t const* _S2,
                                      rsize_t        _N
            );

    

    __inline int __cdecl fwide(
          FILE* _F,
              int   _M
        )
    {
        (void)_F;
        return (_M);
    }

    __inline int __cdecl mbsinit(
          mbstate_t const* _P
        )
    {
        return _P == 0 || _P->_Wchar == 0;
    }

    __inline wchar_t const* __cdecl wmemchr(
          wchar_t const* _S,
                    wchar_t        _C,
                    size_t         _N
        )
    {
        for (; 0 < _N; ++_S, --_N)
            if (*_S == _C)
                return (wchar_t const*)_S;

        return 0;
    }

    __inline int __cdecl wmemcmp(
          wchar_t const* _S1,
          wchar_t const* _S2,
                    size_t         _N
        )
    {
        for (; 0 < _N; ++_S1, ++_S2, --_N)
            if (*_S1 != *_S2)
                return *_S1 < *_S2 ? -1 : 1;

        return 0;
    }

     
    
    __inline 
    wchar_t* __cdecl wmemcpy(
          wchar_t*       _S1,
                wchar_t const* _S2,
                          size_t         _N
        )
    {
        #pragma warning(push)
        #pragma warning(disable : 4995 4996 6386)
        return (wchar_t*)memcpy(_S1, _S2, _N*sizeof(wchar_t));
        #pragma warning(pop)
    }

    __inline 
    wchar_t* __cdecl wmemmove(
          wchar_t*       _S1,
                wchar_t const* _S2,
                              size_t         _N
        )
    {
        #pragma warning(push)
        #pragma warning(disable : 4996 6386)
        return (wchar_t*)memmove(_S1, _S2, _N*sizeof(wchar_t));
        #pragma warning(pop)
    }

     
    
    __inline wchar_t* __cdecl wmemset(
          wchar_t* _S,
                          wchar_t  _C,
                          size_t   _N
        )
    {
        wchar_t *_Su = _S;
        for (; 0 < _N; ++_Su, --_N)
        {
            *_Su = _C;
        }
        return _S;
    }

    

        extern "C++" inline wchar_t* __cdecl wmemchr(
              wchar_t* _S,
                        wchar_t  _C,
                        size_t   _N
            )
        {
            wchar_t const* const _SC = _S;
            return const_cast<wchar_t*>(wmemchr(_SC, _C, _N));
        }

    





} __pragma(pack(pop))



typedef mbstate_t _Mbstatet;

 
namespace std {
using :: _Mbstatet;

using :: mbstate_t; using :: size_t; using :: tm; using :: wint_t;

using :: btowc; using :: fgetwc; using :: fgetws; using :: fputwc;
using :: fputws; using :: fwide; using :: fwprintf;
using :: fwscanf; using :: getwc; using :: getwchar;
using :: mbrlen; using :: mbrtowc; using :: mbsrtowcs;
using :: mbsinit; using :: putwc; using :: putwchar;
using :: swprintf; using :: swscanf; using :: ungetwc;
using :: vfwprintf; using :: vswprintf; using :: vwprintf;
using :: wcrtomb; using :: wprintf; using :: wscanf;
using :: wcsrtombs; using :: wcstol; using :: wcscat;
using :: wcschr; using :: wcscmp; using :: wcscoll;
using :: wcscpy; using :: wcscspn; using :: wcslen;
using :: wcsncat; using :: wcsncmp; using :: wcsncpy;
using :: wcspbrk; using :: wcsrchr; using :: wcsspn;
using :: wcstod; using :: wcstoul; using :: wcsstr;
using :: wcstok; using :: wcsxfrm; using :: wctob;
using :: wmemchr; using :: wmemcmp; using :: wmemcpy;
using :: wmemmove; using :: wmemset; using :: wcsftime;

using :: vfwscanf; using :: vswscanf; using :: vwscanf;
using :: wcstof; using :: wcstold;
using :: wcstoll; using :: wcstoull;
}
 










#pragma once





#pragma once










 


 
namespace std {
using :: ptrdiff_t; using :: size_t;
}
 

 
namespace std {
typedef double max_align_t;	
}

using ::std:: max_align_t;	
 











#pragma once





 #pragma pack(push,8)
 #pragma warning(push,3)
 
 

namespace std {
		
template<class _Elem>
	class initializer_list
	{	
public:
	typedef _Elem value_type;
	typedef const _Elem& reference;
	typedef const _Elem& const_reference;
	typedef size_t size_type;

	typedef const _Elem* iterator;
	typedef const _Elem* const_iterator;

	constexpr initializer_list() noexcept
		: _First(0), _Last(0)
		{	
		}

	constexpr initializer_list(const _Elem *_First_arg,
		const _Elem *_Last_arg) noexcept
		: _First(_First_arg), _Last(_Last_arg)
		{	
		}

	constexpr const _Elem *begin() const noexcept
		{	
		return (_First);
		}

	constexpr const _Elem *end() const noexcept
		{	
		return (_Last);
		}

	constexpr size_t size() const noexcept
		{	
		return ((size_t)(_Last - _First));
		}

private:
	const _Elem *_First;
	const _Elem *_Last;
	};

		
template<class _Elem> inline
	constexpr const _Elem *begin(initializer_list<_Elem> _Ilist) noexcept
	{	
	return (_Ilist.begin());
	}

		
template<class _Elem> inline
	constexpr const _Elem *end(initializer_list<_Elem> _Ilist) noexcept
	{	
	return (_Ilist.end());
	}
}
 
 #pragma warning(pop)
 #pragma pack(pop)











 #pragma pack(push,8)
 #pragma warning(push,3)
 
 

 
  
  
  
 

namespace std {
		
 
 
 

 
 
 
 
 

 
 

  

  












   
   
  

 





















		


		
 

 



































 
 

		

template<bool,
	class _Ty1,
	class _Ty2>
	struct _If
	{	
	typedef _Ty2 type;
	};

template<class _Ty1,
	class _Ty2>
	struct _If<true, _Ty1, _Ty2>
	{	
	typedef _Ty1 type;
	};

template<class _Ty>
	struct _Always_false
	{	
	static constexpr bool value = false;
	};

		

 
		
template<class _Arg,
	class _Result>
	struct unary_function
	{	
	typedef _Arg argument_type;
	typedef _Result result_type;
	};

		
template<class _Arg1,
	class _Arg2,
	class _Result>
	struct binary_function
	{	
	typedef _Arg1 first_argument_type;
	typedef _Arg2 second_argument_type;
	typedef _Result result_type;
	};
 

		
template<class _Ty = void>
	struct plus
	{	
	typedef _Ty first_argument_type;
	typedef _Ty second_argument_type;
	typedef _Ty result_type;

	constexpr _Ty operator()(const _Ty& _Left, const _Ty& _Right) const
		{	
		return (_Left + _Right);
		}
	};

		
template<class _Ty = void>
	struct minus
	{	
	typedef _Ty first_argument_type;
	typedef _Ty second_argument_type;
	typedef _Ty result_type;

	constexpr _Ty operator()(const _Ty& _Left, const _Ty& _Right) const
		{	
		return (_Left - _Right);
		}
	};

		
template<class _Ty = void>
	struct multiplies
	{	
	typedef _Ty first_argument_type;
	typedef _Ty second_argument_type;
	typedef _Ty result_type;

	constexpr _Ty operator()(const _Ty& _Left, const _Ty& _Right) const
		{	
		return (_Left * _Right);
		}
	};

		
template<class _Ty = void>
	struct equal_to
	{	
	typedef _Ty first_argument_type;
	typedef _Ty second_argument_type;
	typedef bool result_type;

	constexpr bool operator()(const _Ty& _Left, const _Ty& _Right) const
		{	
		return (_Left == _Right);
		}
	};

		
template<class _Ty = void>
	struct less
	{	
	typedef _Ty first_argument_type;
	typedef _Ty second_argument_type;
	typedef bool result_type;

	constexpr bool operator()(const _Ty& _Left, const _Ty& _Right) const
		{	
		return (_Left < _Right);
		}
	};

		
template<>
	struct plus<void>
	{	
	typedef int is_transparent;

	template<class _Ty1,
		class _Ty2>
		constexpr auto operator()(_Ty1&& _Left, _Ty2&& _Right) const
		-> decltype(static_cast<_Ty1&&>(_Left)
			+ static_cast<_Ty2&&>(_Right))
		{	
		return (static_cast<_Ty1&&>(_Left)
			+ static_cast<_Ty2&&>(_Right));
		}
	};

		
template<>
	struct minus<void>
	{	
	typedef int is_transparent;

	template<class _Ty1,
		class _Ty2>
		constexpr auto operator()(_Ty1&& _Left, _Ty2&& _Right) const
		-> decltype(static_cast<_Ty1&&>(_Left)
			- static_cast<_Ty2&&>(_Right))
		{	
		return (static_cast<_Ty1&&>(_Left)
			- static_cast<_Ty2&&>(_Right));
		}
	};

		
template<>
	struct multiplies<void>
	{	
	typedef int is_transparent;

	template<class _Ty1,
		class _Ty2>
		constexpr auto operator()(_Ty1&& _Left, _Ty2&& _Right) const
		-> decltype(static_cast<_Ty1&&>(_Left)
			* static_cast<_Ty2&&>(_Right))
		{	
		return (static_cast<_Ty1&&>(_Left)
			* static_cast<_Ty2&&>(_Right));
		}
	};

		
template<>
	struct equal_to<void>
	{	
	typedef int is_transparent;

	template<class _Ty1,
		class _Ty2>
		constexpr auto operator()(_Ty1&& _Left, _Ty2&& _Right) const
		-> decltype(static_cast<_Ty1&&>(_Left)
			== static_cast<_Ty2&&>(_Right))
		{	
		return (static_cast<_Ty1&&>(_Left)
			== static_cast<_Ty2&&>(_Right));
		}
	};

		
template<>
	struct less<void>
	{	
	typedef int is_transparent;

	template<class _Ty1,
		class _Ty2>
		constexpr auto operator()(_Ty1&& _Left, _Ty2&& _Right) const
		-> decltype(static_cast<_Ty1&&>(_Left)
			< static_cast<_Ty2&&>(_Right))
		{	
		return (static_cast<_Ty1&&>(_Left)
			< static_cast<_Ty2&&>(_Right));
		}
	};


}



namespace std {
	
inline size_t _Hash_seq(const unsigned char *_First, size_t _Count)
	{	
 
	static_assert(sizeof(size_t) == 8, "This code is for 64-bit size_t.");
	const size_t _FNV_offset_basis = 14695981039346656037ULL;
	const size_t _FNV_prime = 1099511628211ULL;

 





	size_t _Val = _FNV_offset_basis;
	for (size_t _Next = 0; _Next < _Count; ++_Next)
		{	
		_Val ^= (size_t)_First[_Next];
		_Val *= _FNV_prime;
		}
	return (_Val);
	}

	
template<class _Kty>
	struct _Bitwise_hash
	{	
	typedef _Kty argument_type;
	typedef size_t result_type;

	size_t operator()(const _Kty& _Keyval) const
		{	
		return (_Hash_seq((const unsigned char *)&_Keyval, sizeof (_Kty)));
		}
	};

	
template<class _Kty>
	struct hash
		: public _Bitwise_hash<_Kty>
	{	
	static constexpr bool _Value = __is_enum(_Kty);
	static_assert(_Value,
		"The C++ Standard doesn't provide a hash for this type.");
	};
template<>
	struct hash<bool>
		: public _Bitwise_hash<bool>
	{	
	};

template<>
	struct hash<char>
		: public _Bitwise_hash<char>
	{	
	};

template<>
	struct hash<signed char>
		: public _Bitwise_hash<signed char>
	{	
	};

template<>
	struct hash<unsigned char>
		: public _Bitwise_hash<unsigned char>
	{	
	};

template<>
	struct hash<char16_t>
		: public _Bitwise_hash<char16_t>
	{	
	};

template<>
	struct hash<char32_t>
		: public _Bitwise_hash<char32_t>
	{	
	};

 
template<>
	struct hash<wchar_t>
		: public _Bitwise_hash<wchar_t>
	{	
	};
 

template<>
	struct hash<short>
		: public _Bitwise_hash<short>
	{	
	};

template<>
	struct hash<unsigned short>
		: public _Bitwise_hash<unsigned short>
	{	
	};

template<>
	struct hash<int>
		: public _Bitwise_hash<int>
	{	
	};

template<>
	struct hash<unsigned int>
		: public _Bitwise_hash<unsigned int>
	{	
	};

template<>
	struct hash<long>
		: public _Bitwise_hash<long>
	{	
	};

template<>
	struct hash<unsigned long>
		: public _Bitwise_hash<unsigned long>
	{	
	};

template<>
	struct hash<long long>
		: public _Bitwise_hash<long long>
	{	
	};

template<>
	struct hash<unsigned long long>
		: public _Bitwise_hash<unsigned long long>
	{	
	};

template<>
	struct hash<float>
		: public _Bitwise_hash<float>
	{	
	typedef float _Kty;
	typedef _Bitwise_hash<_Kty> _Mybase;

	size_t operator()(const _Kty& _Keyval) const
		{	
		return (_Mybase::operator()(
			_Keyval == 0 ? 0 : _Keyval));	
		}
	};

template<>
	struct hash<double>
		: public _Bitwise_hash<double>
	{	
	typedef double _Kty;
	typedef _Bitwise_hash<_Kty> _Mybase;

	size_t operator()(const _Kty& _Keyval) const
		{	
		return (_Mybase::operator()(
			_Keyval == 0 ? 0 : _Keyval));	
		}
	};

template<>
	struct hash<long double>
		: public _Bitwise_hash<long double>
	{	
	typedef long double _Kty;
	typedef _Bitwise_hash<_Kty> _Mybase;

	size_t operator()(const _Kty& _Keyval) const
		{	
		return (_Mybase::operator()(
			_Keyval == 0 ? 0 : _Keyval));	
		}
	};

template<class _Ty>
	struct hash<_Ty *>
		: public _Bitwise_hash<_Ty *>
	{	
	};
}



namespace std {
namespace tr1 {	
using ::std:: hash;
}	
}





 





 

 





 

 








 

 



 




































































namespace std {
	
template<class... _Types>
	struct _Arg_types
	{	
	};

template<class _Ty1>
	struct _Arg_types<_Ty1>
	{	
	typedef _Ty1 argument_type;
	};

template<class _Ty1,
	class _Ty2>
	struct _Arg_types<_Ty1, _Ty2>
	{	
	typedef _Ty1 first_argument_type;
	typedef _Ty2 second_argument_type;
	};

	
template<class _Ty>
	struct _Is_function
	{	
	typedef false_type _Bool_type;
	static constexpr bool _Weird = false;
	};












template<class _Ret, class... _Types> struct _Is_function<_Ret __cdecl (_Types...) > : _Arg_types<_Types...> { typedef true_type _Bool_type; typedef _Ret result_type; static constexpr bool _Weird = false; };    template<class _Ret, class... _Types> struct _Is_function<_Ret __vectorcall (_Types...) > : _Arg_types<_Types...> { typedef true_type _Bool_type; typedef _Ret result_type; static constexpr bool _Weird = false; }; template<class _Ret, class... _Types> struct _Is_function<_Ret __cdecl (_Types...) const> : _Arg_types<_Types...> { typedef true_type _Bool_type; typedef _Ret result_type; static constexpr bool _Weird = true; };    template<class _Ret, class... _Types> struct _Is_function<_Ret __vectorcall (_Types...) const> : _Arg_types<_Types...> { typedef true_type _Bool_type; typedef _Ret result_type; static constexpr bool _Weird = true; }; template<class _Ret, class... _Types> struct _Is_function<_Ret __cdecl (_Types...) volatile> : _Arg_types<_Types...> { typedef true_type _Bool_type; typedef _Ret result_type; static constexpr bool _Weird = true; };    template<class _Ret, class... _Types> struct _Is_function<_Ret __vectorcall (_Types...) volatile> : _Arg_types<_Types...> { typedef true_type _Bool_type; typedef _Ret result_type; static constexpr bool _Weird = true; }; template<class _Ret, class... _Types> struct _Is_function<_Ret __cdecl (_Types...) const volatile> : _Arg_types<_Types...> { typedef true_type _Bool_type; typedef _Ret result_type; static constexpr bool _Weird = true; };    template<class _Ret, class... _Types> struct _Is_function<_Ret __vectorcall (_Types...) const volatile> : _Arg_types<_Types...> { typedef true_type _Bool_type; typedef _Ret result_type; static constexpr bool _Weird = true; }; template<class _Ret, class... _Types> struct _Is_function<_Ret __cdecl (_Types...) &> : _Arg_types<_Types...> { typedef true_type _Bool_type; typedef _Ret result_type; static constexpr bool _Weird = true; };    template<class _Ret, class... _Types> struct _Is_function<_Ret __vectorcall (_Types...) &> : _Arg_types<_Types...> { typedef true_type _Bool_type; typedef _Ret result_type; static constexpr bool _Weird = true; }; template<class _Ret, class... _Types> struct _Is_function<_Ret __cdecl (_Types...) const &> : _Arg_types<_Types...> { typedef true_type _Bool_type; typedef _Ret result_type; static constexpr bool _Weird = true; };    template<class _Ret, class... _Types> struct _Is_function<_Ret __vectorcall (_Types...) const &> : _Arg_types<_Types...> { typedef true_type _Bool_type; typedef _Ret result_type; static constexpr bool _Weird = true; }; template<class _Ret, class... _Types> struct _Is_function<_Ret __cdecl (_Types...) volatile &> : _Arg_types<_Types...> { typedef true_type _Bool_type; typedef _Ret result_type; static constexpr bool _Weird = true; };    template<class _Ret, class... _Types> struct _Is_function<_Ret __vectorcall (_Types...) volatile &> : _Arg_types<_Types...> { typedef true_type _Bool_type; typedef _Ret result_type; static constexpr bool _Weird = true; }; template<class _Ret, class... _Types> struct _Is_function<_Ret __cdecl (_Types...) const volatile &> : _Arg_types<_Types...> { typedef true_type _Bool_type; typedef _Ret result_type; static constexpr bool _Weird = true; };    template<class _Ret, class... _Types> struct _Is_function<_Ret __vectorcall (_Types...) const volatile &> : _Arg_types<_Types...> { typedef true_type _Bool_type; typedef _Ret result_type; static constexpr bool _Weird = true; }; template<class _Ret, class... _Types> struct _Is_function<_Ret __cdecl (_Types...) &&> : _Arg_types<_Types...> { typedef true_type _Bool_type; typedef _Ret result_type; static constexpr bool _Weird = true; };    template<class _Ret, class... _Types> struct _Is_function<_Ret __vectorcall (_Types...) &&> : _Arg_types<_Types...> { typedef true_type _Bool_type; typedef _Ret result_type; static constexpr bool _Weird = true; }; template<class _Ret, class... _Types> struct _Is_function<_Ret __cdecl (_Types...) const &&> : _Arg_types<_Types...> { typedef true_type _Bool_type; typedef _Ret result_type; static constexpr bool _Weird = true; };    template<class _Ret, class... _Types> struct _Is_function<_Ret __vectorcall (_Types...) const &&> : _Arg_types<_Types...> { typedef true_type _Bool_type; typedef _Ret result_type; static constexpr bool _Weird = true; }; template<class _Ret, class... _Types> struct _Is_function<_Ret __cdecl (_Types...) volatile &&> : _Arg_types<_Types...> { typedef true_type _Bool_type; typedef _Ret result_type; static constexpr bool _Weird = true; };    template<class _Ret, class... _Types> struct _Is_function<_Ret __vectorcall (_Types...) volatile &&> : _Arg_types<_Types...> { typedef true_type _Bool_type; typedef _Ret result_type; static constexpr bool _Weird = true; }; template<class _Ret, class... _Types> struct _Is_function<_Ret __cdecl (_Types...) const volatile &&> : _Arg_types<_Types...> { typedef true_type _Bool_type; typedef _Ret result_type; static constexpr bool _Weird = true; };    template<class _Ret, class... _Types> struct _Is_function<_Ret __vectorcall (_Types...) const volatile &&> : _Arg_types<_Types...> { typedef true_type _Bool_type; typedef _Ret result_type; static constexpr bool _Weird = true; };












template<class _Ret, class... _Types> struct _Is_function<_Ret (_Types..., ...) > { typedef true_type _Bool_type; typedef _Ret result_type; static constexpr bool _Weird = false; }; template<class _Ret, class... _Types> struct _Is_function<_Ret (_Types..., ...) const> { typedef true_type _Bool_type; typedef _Ret result_type; static constexpr bool _Weird = true; }; template<class _Ret, class... _Types> struct _Is_function<_Ret (_Types..., ...) volatile> { typedef true_type _Bool_type; typedef _Ret result_type; static constexpr bool _Weird = true; }; template<class _Ret, class... _Types> struct _Is_function<_Ret (_Types..., ...) const volatile> { typedef true_type _Bool_type; typedef _Ret result_type; static constexpr bool _Weird = true; }; template<class _Ret, class... _Types> struct _Is_function<_Ret (_Types..., ...) &> { typedef true_type _Bool_type; typedef _Ret result_type; static constexpr bool _Weird = true; }; template<class _Ret, class... _Types> struct _Is_function<_Ret (_Types..., ...) const &> { typedef true_type _Bool_type; typedef _Ret result_type; static constexpr bool _Weird = true; }; template<class _Ret, class... _Types> struct _Is_function<_Ret (_Types..., ...) volatile &> { typedef true_type _Bool_type; typedef _Ret result_type; static constexpr bool _Weird = true; }; template<class _Ret, class... _Types> struct _Is_function<_Ret (_Types..., ...) const volatile &> { typedef true_type _Bool_type; typedef _Ret result_type; static constexpr bool _Weird = true; }; template<class _Ret, class... _Types> struct _Is_function<_Ret (_Types..., ...) &&> { typedef true_type _Bool_type; typedef _Ret result_type; static constexpr bool _Weird = true; }; template<class _Ret, class... _Types> struct _Is_function<_Ret (_Types..., ...) const &&> { typedef true_type _Bool_type; typedef _Ret result_type; static constexpr bool _Weird = true; }; template<class _Ret, class... _Types> struct _Is_function<_Ret (_Types..., ...) volatile &&> { typedef true_type _Bool_type; typedef _Ret result_type; static constexpr bool _Weird = true; }; template<class _Ret, class... _Types> struct _Is_function<_Ret (_Types..., ...) const volatile &&> { typedef true_type _Bool_type; typedef _Ret result_type; static constexpr bool _Weird = true; };


template<class _Ty>
	struct is_function
		: _Is_function<_Ty>::_Bool_type
	{	
	};

 
template<class _Ty>
	constexpr bool is_function_v = is_function<_Ty>::value;
 

		





















template<class _Ty> inline
	constexpr _Ty *addressof(_Ty& _Val) noexcept
	{	
	return (__builtin_addressof(_Val));
	}



		
template<class _Ptrty> inline
	auto _Unfancy(_Ptrty _Ptr)
	{	
	return (::std:: addressof(*_Ptr));
	}

template<class _Ty> inline
	_Ty * _Unfancy(_Ty * _Ptr)
	{	
	return (_Ptr);
	}

}
 
 #pragma warning(pop)
 #pragma pack(pop)









 #pragma pack(push,8)
 #pragma warning(push,3)
 
 

namespace std {















 

 
 
 
 

 
 
 
 

 
 
 
 

		
typedef enum
	{	
	denorm_indeterminate = -1,
	denorm_absent = 0,
	denorm_present = 1}
		float_denorm_style;

		
typedef enum
	{	
	round_indeterminate = -1,
	round_toward_zero = 0,
	round_to_nearest = 1,
	round_toward_infinity = 2,
	round_toward_neg_infinity = 3}
		float_round_style;

		
struct _Num_base
	{	
	static constexpr float_denorm_style has_denorm = (float_denorm_style)(denorm_absent);
	static constexpr bool has_denorm_loss = (bool)(false);
	static constexpr bool has_infinity = (bool)(false);
	static constexpr bool has_quiet_NaN = (bool)(false);
	static constexpr bool has_signaling_NaN = (bool)(false);
	static constexpr bool is_bounded = (bool)(false);
	static constexpr bool is_exact = (bool)(false);
	static constexpr bool is_iec559 = (bool)(false);
	static constexpr bool is_integer = (bool)(false);
	static constexpr bool is_modulo = (bool)(false);
	static constexpr bool is_signed = (bool)(false);
	static constexpr bool is_specialized = (bool)(false);
	static constexpr bool tinyness_before = (bool)(false);
	static constexpr bool traps = (bool)(false);
	static constexpr float_round_style round_style = (float_round_style)(round_toward_zero);
	static constexpr int digits = (int)(0);
	static constexpr int digits10 = (int)(0);

	static constexpr int max_digits10 = (int)(0);

	static constexpr int max_exponent = (int)(0);
	static constexpr int max_exponent10 = (int)(0);
	static constexpr int min_exponent = (int)(0);
	static constexpr int min_exponent10 = (int)(0);
	static constexpr int radix = (int)(0);
	};

		
template<class _Ty>
	class numeric_limits
		: public _Num_base
	{	
public:
	static constexpr _Ty (min)() noexcept
		{	
		return (_Ty());
		}

	static constexpr _Ty (max)() noexcept
		{	
		return (_Ty());
		}

	static constexpr _Ty lowest() noexcept
		{	
		return (_Ty());
		}

	static constexpr _Ty epsilon() noexcept
		{	
		return (_Ty());
		}

	static constexpr _Ty round_error() noexcept
		{	
		return (_Ty());
		}

	static constexpr _Ty denorm_min() noexcept
		{	
		return (_Ty());
		}

	static constexpr _Ty infinity() noexcept
		{	
		return (_Ty());
		}

	static constexpr _Ty quiet_NaN() noexcept
		{	
		return (_Ty());
		}

	static constexpr _Ty signaling_NaN() noexcept
		{	
		return (_Ty());
		}
	};

template<class _Ty>
	class numeric_limits<const _Ty>
		: public numeric_limits<_Ty>
	{	
	};

template<class _Ty>
	class numeric_limits<volatile _Ty>
		: public numeric_limits<_Ty>
	{	
	};

template<class _Ty>
	class numeric_limits<const volatile _Ty>
		: public numeric_limits<_Ty>
	{	
	};

		
struct _Num_int_base
	: public _Num_base
	{	
	static constexpr bool is_bounded = (bool)(true);
	static constexpr bool is_exact = (bool)(true);
	static constexpr bool is_integer = (bool)(true);
	static constexpr bool is_modulo = (bool)(true);
	static constexpr bool is_specialized = (bool)(true);
	static constexpr int radix = (int)(2);
	};

		
struct _Num_float_base
	: public _Num_base
	{	
	static constexpr float_denorm_style has_denorm = (float_denorm_style)(denorm_present);
	static constexpr bool has_denorm_loss = (bool)(true);
	static constexpr bool has_infinity = (bool)(true);
	static constexpr bool has_quiet_NaN = (bool)(true);
	static constexpr bool has_signaling_NaN = (bool)(true);
	static constexpr bool is_bounded = (bool)(true);
	static constexpr bool is_exact = (bool)(false);
	static constexpr bool is_iec559 = (bool)(true);
	static constexpr bool is_integer = (bool)(false);
	static constexpr bool is_modulo = (bool)(false);
	static constexpr bool is_signed = (bool)(true);
	static constexpr bool is_specialized = (bool)(true);
	static constexpr bool tinyness_before = (bool)(true);
	static constexpr bool traps = (bool)(false);
	static constexpr float_round_style round_style = (float_round_style)(round_to_nearest);
	static constexpr int radix = (int)(2);
	};

		
template<> class numeric_limits<char>
	: public _Num_int_base
	{	
public:
	typedef char _Ty;

	static constexpr _Ty (min)() noexcept
		{	
		return ((-128));
		}

	static constexpr _Ty (max)() noexcept
		{	
		return (127);
		}

	static constexpr _Ty lowest() noexcept
		{	
		return ((min)());
		}

	static constexpr _Ty epsilon() noexcept
		{	
		return (0);
		}

	static constexpr _Ty round_error() noexcept
		{	
		return (0);
		}

	static constexpr _Ty denorm_min() noexcept
		{	
		return (0);
		}

	static constexpr _Ty infinity() noexcept
		{	
		return (0);
		}

	static constexpr _Ty quiet_NaN() noexcept
		{	
		return (0);
		}

	static constexpr _Ty signaling_NaN() noexcept
		{	
		return (0);
		}

	static constexpr bool is_signed = (bool)((-128) != 0);
	static constexpr int digits = (int)(8 - ((-128) != 0 ? 1 : 0));
	static constexpr int digits10 = (int)((8 - ((-128) != 0 ? 1 : 0)) * 301L / 1000);
	};

		
template<> class numeric_limits<wchar_t>
	: public _Num_int_base
	{	
public:
	typedef wchar_t _Ty;

	static constexpr _Ty (min)() noexcept
		{	
		return ((_Ty)0x0000);
		}

	static constexpr _Ty (max)() noexcept
		{	
		return ((_Ty)0xffff);
		}

	static constexpr _Ty lowest() noexcept
		{	
		return ((min)());
		}

	static constexpr _Ty epsilon() noexcept
		{	
		return (0);
		}

	static constexpr _Ty round_error() noexcept
		{	
		return (0);
		}

	static constexpr _Ty denorm_min() noexcept
		{	
		return (0);
		}

	static constexpr _Ty infinity() noexcept
		{	
		return (0);
		}

	static constexpr _Ty quiet_NaN() noexcept
		{	
		return (0);
		}

	static constexpr _Ty signaling_NaN() noexcept
		{	
		return (0);
		}

	static constexpr bool is_signed = (bool)(0x0000 != 0);
	static constexpr int digits = (int)(8 * sizeof (wchar_t) - (0x0000 != 0 ? 1 : 0));
	static constexpr int digits10 = (int)((8 * sizeof (wchar_t) - (0x0000 != 0 ? 1 : 0)) * 301L / 1000);
	};

		
template<> class numeric_limits<bool>
	: public _Num_int_base
	{	
public:
	typedef bool _Ty;

	static constexpr _Ty (min)() noexcept
		{	
		return (false);
		}

	static constexpr _Ty (max)() noexcept
		{	
		return (true);
		}

	static constexpr _Ty lowest() noexcept
		{	
		return ((min)());
		}

	static constexpr _Ty epsilon() noexcept
		{	
		return (0);
		}

	static constexpr _Ty round_error() noexcept
		{	
		return (0);
		}

	static constexpr _Ty denorm_min() noexcept
		{	
		return (0);
		}

	static constexpr _Ty infinity() noexcept
		{	
		return (0);
		}

	static constexpr _Ty quiet_NaN() noexcept
		{	
		return (0);
		}

	static constexpr _Ty signaling_NaN() noexcept
		{	
		return (0);
		}

	static constexpr bool is_modulo = (bool)(false);
	static constexpr bool is_signed = (bool)(false);
	static constexpr int digits = (int)(1);
	static constexpr int digits10 = (int)(0);
	};

		
template<> class numeric_limits<signed char>
	: public _Num_int_base
	{	
public:
	typedef signed char _Ty;

	static constexpr _Ty (min)() noexcept
		{	
		return ((-128));
		}

	static constexpr _Ty (max)() noexcept
		{	
		return (127);
		}

	static constexpr _Ty lowest() noexcept
		{	
		return ((min)());
		}

	static constexpr _Ty epsilon() noexcept
		{	
		return (0);
		}

	static constexpr _Ty round_error() noexcept
		{	
		return (0);
		}

	static constexpr _Ty denorm_min() noexcept
		{	
		return (0);
		}

	static constexpr _Ty infinity() noexcept
		{	
		return (0);
		}

	static constexpr _Ty quiet_NaN() noexcept
		{	
		return (0);
		}

	static constexpr _Ty signaling_NaN() noexcept
		{	
		return (0);
		}

	static constexpr bool is_signed = (bool)(true);
	static constexpr int digits = (int)(8 - 1);
	static constexpr int digits10 = (int)((8 - 1) * 301L / 1000);
	};

		
template<> class numeric_limits<unsigned char>
	: public _Num_int_base
	{	
public:
	typedef unsigned char _Ty;

	static constexpr _Ty (min)() noexcept
		{	
		return (0);
		}

	static constexpr _Ty (max)() noexcept
		{	
		return (0xff);
		}

	static constexpr _Ty lowest() noexcept
		{	
		return ((min)());
		}

	static constexpr _Ty epsilon() noexcept
		{	
		return (0);
		}

	static constexpr _Ty round_error() noexcept
		{	
		return (0);
		}

	static constexpr _Ty denorm_min() noexcept
		{	
		return (0);
		}

	static constexpr _Ty infinity() noexcept
		{	
		return (0);
		}

	static constexpr _Ty quiet_NaN() noexcept
		{	
		return (0);
		}

	static constexpr _Ty signaling_NaN() noexcept
		{	
		return (0);
		}

	static constexpr bool is_signed = (bool)(false);
	static constexpr int digits = (int)(8);
	static constexpr int digits10 = (int)(8 * 301L / 1000);
	};

		
template<> class numeric_limits<short>
	: public _Num_int_base
	{	
public:
	typedef short _Ty;

	static constexpr _Ty (min)() noexcept
		{	
		return ((-32768));
		}

	static constexpr _Ty (max)() noexcept
		{	
		return (32767);
		}

	static constexpr _Ty lowest() noexcept
		{	
		return ((min)());
		}

	static constexpr _Ty epsilon() noexcept
		{	
		return (0);
		}

	static constexpr _Ty round_error() noexcept
		{	
		return (0);
		}

	static constexpr _Ty denorm_min() noexcept
		{	
		return (0);
		}

	static constexpr _Ty infinity() noexcept
		{	
		return (0);
		}

	static constexpr _Ty quiet_NaN() noexcept
		{	
		return (0);
		}

	static constexpr _Ty signaling_NaN() noexcept
		{	
		return (0);
		}

	static constexpr bool is_signed = (bool)(true);
	static constexpr int digits = (int)(8 * sizeof (short) - 1);
	static constexpr int digits10 = (int)((8 * sizeof (short) - 1) * 301L / 1000);
	};

 
		
template<> class numeric_limits<unsigned short>
	: public _Num_int_base
	{	
public:
	typedef unsigned short _Ty;

	static constexpr _Ty (min)() noexcept
		{	
		return (0);
		}

	static constexpr _Ty (max)() noexcept
		{	
		return (0xffff);
		}

	static constexpr _Ty lowest() noexcept
		{	
		return ((min)());
		}

	static constexpr _Ty epsilon() noexcept
		{	
		return (0);
		}

	static constexpr _Ty round_error() noexcept
		{	
		return (0);
		}

	static constexpr _Ty denorm_min() noexcept
		{	
		return (0);
		}

	static constexpr _Ty infinity() noexcept
		{	
		return (0);
		}

	static constexpr _Ty quiet_NaN() noexcept
		{	
		return (0);
		}

	static constexpr _Ty signaling_NaN() noexcept
		{	
		return (0);
		}

	static constexpr bool is_signed = (bool)(false);
	static constexpr int digits = (int)(8 * sizeof (unsigned short));
	static constexpr int digits10 = (int)(8 * sizeof (unsigned short) * 301L / 1000);
	};
 

		
template<> class numeric_limits<char16_t>
	: public _Num_int_base
	{	
public:
	typedef char16_t _Ty;

	static constexpr _Ty (min)() noexcept
		{	
		return (0);
		}

	static constexpr _Ty (max)() noexcept
		{	
		return (0xffff);
		}

	static constexpr _Ty lowest() noexcept
		{	
		return ((min)());
		}

	static constexpr _Ty epsilon() noexcept
		{	
		return (0);
		}

	static constexpr _Ty round_error() noexcept
		{	
		return (0);
		}

	static constexpr _Ty denorm_min() noexcept
		{	
		return (0);
		}

	static constexpr _Ty infinity() noexcept
		{	
		return (0);
		}

	static constexpr _Ty quiet_NaN() noexcept
		{	
		return (0);
		}

	static constexpr _Ty signaling_NaN() noexcept
		{	
		return (0);
		}

	static constexpr bool is_signed = (bool)(false);
	static constexpr int digits = (int)(8 * sizeof (char16_t));
	static constexpr int digits10 = (int)(8 * sizeof (char16_t) * 301L / 1000);
	};

		
template<> class numeric_limits<int>
	: public _Num_int_base
	{	
public:
	typedef int _Ty;

	static constexpr _Ty (min)() noexcept
		{	
		return ((-2147483647 - 1));
		}

	static constexpr _Ty (max)() noexcept
		{	
		return (2147483647);
		}

	static constexpr _Ty lowest() noexcept
		{	
		return ((min)());
		}

	static constexpr _Ty epsilon() noexcept
		{	
		return (0);
		}

	static constexpr _Ty round_error() noexcept
		{	
		return (0);
		}

	static constexpr _Ty denorm_min() noexcept
		{	
		return (0);
		}

	static constexpr _Ty infinity() noexcept
		{	
		return (0);
		}

	static constexpr _Ty quiet_NaN() noexcept
		{	
		return (0);
		}

	static constexpr _Ty signaling_NaN() noexcept
		{	
		return (0);
		}

	static constexpr bool is_signed = (bool)(true);
	static constexpr int digits = (int)(8 * sizeof (int) - 1);
	static constexpr int digits10 = (int)((8 * sizeof (int) - 1) * 301L / 1000);
	};

		
template<> class numeric_limits<unsigned int>
	: public _Num_int_base
	{	
public:
	typedef unsigned int _Ty;

	static constexpr _Ty (min)() noexcept
		{	
		return (0);
		}

	static constexpr _Ty (max)() noexcept
		{	
		return (0xffffffff);
		}

	static constexpr _Ty lowest() noexcept
		{	
		return ((min)());
		}

	static constexpr _Ty epsilon() noexcept
		{	
		return (0);
		}

	static constexpr _Ty round_error() noexcept
		{	
		return (0);
		}

	static constexpr _Ty denorm_min() noexcept
		{	
		return (0);
		}

	static constexpr _Ty infinity() noexcept
		{	
		return (0);
		}

	static constexpr _Ty quiet_NaN() noexcept
		{	
		return (0);
		}

	static constexpr _Ty signaling_NaN() noexcept
		{	
		return (0);
		}

	static constexpr bool is_signed = (bool)(false);
	static constexpr int digits = (int)(8 * sizeof (unsigned int));
	static constexpr int digits10 = (int)(8 * sizeof (unsigned int) * 301L / 1000);
	};

		
template<> class numeric_limits<long>
	: public _Num_int_base
	{	
public:
	typedef long _Ty;

	static constexpr _Ty (min)() noexcept
		{	
		return ((-2147483647L - 1));
		}

	static constexpr _Ty (max)() noexcept
		{	
		return (2147483647L);
		}

	static constexpr _Ty lowest() noexcept
		{	
		return ((min)());
		}

	static constexpr _Ty epsilon() noexcept
		{	
		return (0);
		}

	static constexpr _Ty round_error() noexcept
		{	
		return (0);
		}

	static constexpr _Ty denorm_min() noexcept
		{	
		return (0);
		}

	static constexpr _Ty infinity() noexcept
		{	
		return (0);
		}

	static constexpr _Ty quiet_NaN() noexcept
		{	
		return (0);
		}

	static constexpr _Ty signaling_NaN() noexcept
		{	
		return (0);
		}

	static constexpr bool is_signed = (bool)(true);
	static constexpr int digits = (int)(8 * sizeof (long) - 1);
	static constexpr int digits10 = (int)((8 * sizeof (long) - 1) * 301L / 1000);
	};

		
template<> class numeric_limits<unsigned long>
	: public _Num_int_base
	{	
public:
	typedef unsigned long _Ty;

	static constexpr _Ty (min)() noexcept
		{	
		return (0);
		}

	static constexpr _Ty (max)() noexcept
		{	
		return (0xffffffffUL);
		}

	static constexpr _Ty lowest() noexcept
		{	
		return ((min)());
		}

	static constexpr _Ty epsilon() noexcept
		{	
		return (0);
		}

	static constexpr _Ty round_error() noexcept
		{	
		return (0);
		}

	static constexpr _Ty denorm_min() noexcept
		{	
		return (0);
		}

	static constexpr _Ty infinity() noexcept
		{	
		return (0);
		}

	static constexpr _Ty quiet_NaN() noexcept
		{	
		return (0);
		}

	static constexpr _Ty signaling_NaN() noexcept
		{	
		return (0);
		}

	static constexpr bool is_signed = (bool)(false);
	static constexpr int digits = (int)(8 * sizeof (unsigned long));
	static constexpr int digits10 = (int)(8 * sizeof (unsigned long) * 301L / 1000);
	};

		
template<> class numeric_limits<char32_t>
	: public _Num_int_base
	{	
public:
	typedef char32_t _Ty;

	static constexpr _Ty (min)() noexcept
		{	
		return (0);
		}

	static constexpr _Ty (max)() noexcept
		{	
		return (0xffffffff);
		}

	static constexpr _Ty lowest() noexcept
		{	
		return ((min)());
		}

	static constexpr _Ty epsilon() noexcept
		{	
		return (0);
		}

	static constexpr _Ty round_error() noexcept
		{	
		return (0);
		}

	static constexpr _Ty denorm_min() noexcept
		{	
		return (0);
		}

	static constexpr _Ty infinity() noexcept
		{	
		return (0);
		}

	static constexpr _Ty quiet_NaN() noexcept
		{	
		return (0);
		}

	static constexpr _Ty signaling_NaN() noexcept
		{	
		return (0);
		}

	static constexpr bool is_signed = (bool)(false);
	static constexpr int digits = (int)(8 * sizeof (char32_t));
	static constexpr int digits10 = (int)(8 * sizeof (char32_t) * 301L / 1000);
	};

		
template<> class numeric_limits<long long>
	: public _Num_int_base
	{	
public:
	typedef long long _Ty;

	static constexpr _Ty (min)() noexcept
		{	
		return (-0x7fffffffffffffff - 1);
		}

	static constexpr _Ty (max)() noexcept
		{	
		return (0x7fffffffffffffff);
		}

	static constexpr _Ty lowest() noexcept
		{	
		return ((min)());
		}

	static constexpr _Ty epsilon() noexcept
		{	
		return (0);
		}

	static constexpr _Ty round_error() noexcept
		{	
		return (0);
		}

	static constexpr _Ty denorm_min() noexcept
		{	
		return (0);
		}

	static constexpr _Ty infinity() noexcept
		{	
		return (0);
		}

	static constexpr _Ty quiet_NaN() noexcept
		{	
		return (0);
		}

	static constexpr _Ty signaling_NaN() noexcept
		{	
		return (0);
		}

	static constexpr bool is_signed = (bool)(true);
	static constexpr int digits = (int)(8 * sizeof (long long) - 1);
	static constexpr int digits10 = (int)((8 * sizeof (long long) - 1) * 301L / 1000);
	};

		
template<> class numeric_limits<unsigned long long>
	: public _Num_int_base
	{	
public:
	typedef unsigned long long _Ty;

	static constexpr _Ty (min)() noexcept
		{	
		return (0);
		}

	static constexpr _Ty (max)() noexcept
		{	
		return (0xffffffffffffffff);
		}

	static constexpr _Ty lowest() noexcept
		{	
		return ((min)());
		}

	static constexpr _Ty epsilon() noexcept
		{	
		return (0);
		}

	static constexpr _Ty round_error() noexcept
		{	
		return (0);
		}

	static constexpr _Ty denorm_min() noexcept
		{	
		return (0);
		}

	static constexpr _Ty infinity() noexcept
		{	
		return (0);
		}

	static constexpr _Ty quiet_NaN() noexcept
		{	
		return (0);
		}

	static constexpr _Ty signaling_NaN() noexcept
		{	
		return (0);
		}

	static constexpr bool is_signed = (bool)(false);
	static constexpr int digits = (int)(8 * sizeof (unsigned long long));
	static constexpr int digits10 = (int)(8 * sizeof (unsigned long long) * 301L / 1000);
	};

		
template<> class numeric_limits<float>
	: public _Num_float_base
	{	
public:
	typedef float _Ty;

	static constexpr _Ty (min)() noexcept
		{	
		return (1.175494351e-38F);
		}

	static constexpr _Ty (max)() noexcept
		{	
		return (3.402823466e+38F);
		}

	static constexpr _Ty lowest() noexcept
		{	
		return (-(max)());
		}

	static constexpr _Ty epsilon() noexcept
		{	
		return (1.192092896e-07F);
		}

	static constexpr _Ty round_error() noexcept
		{	
		return (0.5F);
		}

	static constexpr _Ty denorm_min() noexcept
		{	
		return (1.401298464e-45F);
		}

	static constexpr _Ty infinity() noexcept
		{	
		return (__builtin_huge_valf());
		}

	static constexpr _Ty quiet_NaN() noexcept
		{	
		return (__builtin_nanf("0"));
		}

	static constexpr _Ty signaling_NaN() noexcept
		{	
		return (__builtin_nansf("1"));
		}

	static constexpr int digits = (int)(24);
	static constexpr int digits10 = (int)(6);

	static constexpr int max_digits10 = (int)(2 + 24 * 301L / 1000);

	static constexpr int max_exponent = (int)((int)128);
	static constexpr int max_exponent10 = (int)((int)38);
	static constexpr int min_exponent = (int)((int)(-125));
	static constexpr int min_exponent10 = (int)((int)(-37));
	};

		
template<> class numeric_limits<double>
	: public _Num_float_base
	{	
public:
	typedef double _Ty;

	static constexpr _Ty (min)() noexcept
		{	
		return (2.2250738585072014e-308);
		}

	static constexpr _Ty (max)() noexcept
		{	
		return (1.7976931348623158e+308);
		}

	static constexpr _Ty lowest() noexcept
		{	
		return (-(max)());
		}

	static constexpr _Ty epsilon() noexcept
		{	
		return (2.2204460492503131e-016);
		}

	static constexpr _Ty round_error() noexcept
		{	
		return (0.5);
		}

	static constexpr _Ty denorm_min() noexcept
		{	
		return (4.9406564584124654e-324);
		}

	static constexpr _Ty infinity() noexcept
		{	
		return (__builtin_huge_val());
		}

	static constexpr _Ty quiet_NaN() noexcept
		{	
		return (__builtin_nan("0"));
		}

	static constexpr _Ty signaling_NaN() noexcept
		{	
		return (__builtin_nans("1"));
		}

	static constexpr int digits = (int)(53);
	static constexpr int digits10 = (int)(15);

	static constexpr int max_digits10 = (int)(2 + 53 * 301L / 1000);

	static constexpr int max_exponent = (int)((int)1024);
	static constexpr int max_exponent10 = (int)((int)308);
	static constexpr int min_exponent = (int)((int)(-1021));
	static constexpr int min_exponent10 = (int)((int)(-307));
	};

		
template<> class numeric_limits<long double>
	: public _Num_float_base
	{	
public:
	typedef long double _Ty;

	static constexpr _Ty (min)() noexcept
		{	
		return (2.2250738585072014e-308);
		}

	static constexpr _Ty (max)() noexcept
		{	
		return (1.7976931348623158e+308);
		}

	static constexpr _Ty lowest() noexcept
		{	
		return (-(max)());
		}

	static constexpr _Ty epsilon() noexcept
		{	
		return (2.2204460492503131e-016);
		}

	static constexpr _Ty round_error() noexcept
		{	
		return (0.5L);
		}

	static constexpr _Ty denorm_min() noexcept
		{	
		return (4.9406564584124654e-324);
		}

	static constexpr _Ty infinity() noexcept
		{	
		return (__builtin_huge_val());
		}

	static constexpr _Ty quiet_NaN() noexcept
		{	
		return (__builtin_nan("0"));
		}

	static constexpr _Ty signaling_NaN() noexcept
		{	
		return (__builtin_nans("1"));
		}

	static constexpr int digits = (int)(53);
	static constexpr int digits10 = (int)(15);

	static constexpr int max_digits10 = (int)(2 + 53 * 301L / 1000);

	static constexpr int max_exponent = (int)((int)1024);
	static constexpr int max_exponent10 = (int)((int)308);
	static constexpr int min_exponent = (int)((int)(-1021));
	static constexpr int min_exponent10 = (int)((int)(-307));
	};

  











































































































































































































 
 
 
 

 
 
 
 

 
 
 
 
}
 
 #pragma warning(pop)
 #pragma pack(pop)











#pragma once





#pragma once






#pragma once





 #pragma pack(push,8)
 #pragma warning(push,3)
 
 
 #pragma warning(disable: 4180 4296)

namespace std {
template<class _Ty>
	struct _Is_memfunptr
	{	
	typedef false_type _Bool_type;
	};













template<class _Ret, class _Arg0, class... _Types> struct _Is_memfunptr<_Ret (__cdecl _Arg0::*)(_Types...)  > : _Arg_types< _Arg0 *, _Types...> { typedef true_type _Bool_type; typedef _Ret result_type; typedef _Arg0 _Class_type; };     template<class _Ret, class _Arg0, class... _Types> struct _Is_memfunptr<_Ret (__vectorcall _Arg0::*)(_Types...)  > : _Arg_types< _Arg0 *, _Types...> { typedef true_type _Bool_type; typedef _Ret result_type; typedef _Arg0 _Class_type; }; template<class _Ret, class _Arg0, class... _Types> struct _Is_memfunptr<_Ret (__cdecl _Arg0::*)(_Types...) const > : _Arg_types<const _Arg0 *, _Types...> { typedef true_type _Bool_type; typedef _Ret result_type; typedef _Arg0 _Class_type; };     template<class _Ret, class _Arg0, class... _Types> struct _Is_memfunptr<_Ret (__vectorcall _Arg0::*)(_Types...) const > : _Arg_types<const _Arg0 *, _Types...> { typedef true_type _Bool_type; typedef _Ret result_type; typedef _Arg0 _Class_type; }; template<class _Ret, class _Arg0, class... _Types> struct _Is_memfunptr<_Ret (__cdecl _Arg0::*)(_Types...) volatile > : _Arg_types<volatile _Arg0 *, _Types...> { typedef true_type _Bool_type; typedef _Ret result_type; typedef _Arg0 _Class_type; };     template<class _Ret, class _Arg0, class... _Types> struct _Is_memfunptr<_Ret (__vectorcall _Arg0::*)(_Types...) volatile > : _Arg_types<volatile _Arg0 *, _Types...> { typedef true_type _Bool_type; typedef _Ret result_type; typedef _Arg0 _Class_type; }; template<class _Ret, class _Arg0, class... _Types> struct _Is_memfunptr<_Ret (__cdecl _Arg0::*)(_Types...) const volatile > : _Arg_types<const volatile _Arg0 *, _Types...> { typedef true_type _Bool_type; typedef _Ret result_type; typedef _Arg0 _Class_type; };     template<class _Ret, class _Arg0, class... _Types> struct _Is_memfunptr<_Ret (__vectorcall _Arg0::*)(_Types...) const volatile > : _Arg_types<const volatile _Arg0 *, _Types...> { typedef true_type _Bool_type; typedef _Ret result_type; typedef _Arg0 _Class_type; }; template<class _Ret, class _Arg0, class... _Types> struct _Is_memfunptr<_Ret (__cdecl _Arg0::*)(_Types...)  &> : _Arg_types< _Arg0 *, _Types...> { typedef true_type _Bool_type; typedef _Ret result_type; typedef _Arg0 _Class_type; };     template<class _Ret, class _Arg0, class... _Types> struct _Is_memfunptr<_Ret (__vectorcall _Arg0::*)(_Types...)  &> : _Arg_types< _Arg0 *, _Types...> { typedef true_type _Bool_type; typedef _Ret result_type; typedef _Arg0 _Class_type; }; template<class _Ret, class _Arg0, class... _Types> struct _Is_memfunptr<_Ret (__cdecl _Arg0::*)(_Types...) const &> : _Arg_types<const _Arg0 *, _Types...> { typedef true_type _Bool_type; typedef _Ret result_type; typedef _Arg0 _Class_type; };     template<class _Ret, class _Arg0, class... _Types> struct _Is_memfunptr<_Ret (__vectorcall _Arg0::*)(_Types...) const &> : _Arg_types<const _Arg0 *, _Types...> { typedef true_type _Bool_type; typedef _Ret result_type; typedef _Arg0 _Class_type; }; template<class _Ret, class _Arg0, class... _Types> struct _Is_memfunptr<_Ret (__cdecl _Arg0::*)(_Types...) volatile &> : _Arg_types<volatile _Arg0 *, _Types...> { typedef true_type _Bool_type; typedef _Ret result_type; typedef _Arg0 _Class_type; };     template<class _Ret, class _Arg0, class... _Types> struct _Is_memfunptr<_Ret (__vectorcall _Arg0::*)(_Types...) volatile &> : _Arg_types<volatile _Arg0 *, _Types...> { typedef true_type _Bool_type; typedef _Ret result_type; typedef _Arg0 _Class_type; }; template<class _Ret, class _Arg0, class... _Types> struct _Is_memfunptr<_Ret (__cdecl _Arg0::*)(_Types...) const volatile &> : _Arg_types<const volatile _Arg0 *, _Types...> { typedef true_type _Bool_type; typedef _Ret result_type; typedef _Arg0 _Class_type; };     template<class _Ret, class _Arg0, class... _Types> struct _Is_memfunptr<_Ret (__vectorcall _Arg0::*)(_Types...) const volatile &> : _Arg_types<const volatile _Arg0 *, _Types...> { typedef true_type _Bool_type; typedef _Ret result_type; typedef _Arg0 _Class_type; }; template<class _Ret, class _Arg0, class... _Types> struct _Is_memfunptr<_Ret (__cdecl _Arg0::*)(_Types...)  &&> : _Arg_types< _Arg0 *, _Types...> { typedef true_type _Bool_type; typedef _Ret result_type; typedef _Arg0 _Class_type; };     template<class _Ret, class _Arg0, class... _Types> struct _Is_memfunptr<_Ret (__vectorcall _Arg0::*)(_Types...)  &&> : _Arg_types< _Arg0 *, _Types...> { typedef true_type _Bool_type; typedef _Ret result_type; typedef _Arg0 _Class_type; }; template<class _Ret, class _Arg0, class... _Types> struct _Is_memfunptr<_Ret (__cdecl _Arg0::*)(_Types...) const &&> : _Arg_types<const _Arg0 *, _Types...> { typedef true_type _Bool_type; typedef _Ret result_type; typedef _Arg0 _Class_type; };     template<class _Ret, class _Arg0, class... _Types> struct _Is_memfunptr<_Ret (__vectorcall _Arg0::*)(_Types...) const &&> : _Arg_types<const _Arg0 *, _Types...> { typedef true_type _Bool_type; typedef _Ret result_type; typedef _Arg0 _Class_type; }; template<class _Ret, class _Arg0, class... _Types> struct _Is_memfunptr<_Ret (__cdecl _Arg0::*)(_Types...) volatile &&> : _Arg_types<volatile _Arg0 *, _Types...> { typedef true_type _Bool_type; typedef _Ret result_type; typedef _Arg0 _Class_type; };     template<class _Ret, class _Arg0, class... _Types> struct _Is_memfunptr<_Ret (__vectorcall _Arg0::*)(_Types...) volatile &&> : _Arg_types<volatile _Arg0 *, _Types...> { typedef true_type _Bool_type; typedef _Ret result_type; typedef _Arg0 _Class_type; }; template<class _Ret, class _Arg0, class... _Types> struct _Is_memfunptr<_Ret (__cdecl _Arg0::*)(_Types...) const volatile &&> : _Arg_types<const volatile _Arg0 *, _Types...> { typedef true_type _Bool_type; typedef _Ret result_type; typedef _Arg0 _Class_type; };     template<class _Ret, class _Arg0, class... _Types> struct _Is_memfunptr<_Ret (__vectorcall _Arg0::*)(_Types...) const volatile &&> : _Arg_types<const volatile _Arg0 *, _Types...> { typedef true_type _Bool_type; typedef _Ret result_type; typedef _Arg0 _Class_type; };













template<class _Ret, class _Arg0, class... _Types> struct _Is_memfunptr<_Ret (_Arg0::*)(_Types..., ...) > { typedef true_type _Bool_type; typedef _Ret result_type; typedef _Arg0 _Class_type; }; template<class _Ret, class _Arg0, class... _Types> struct _Is_memfunptr<_Ret (_Arg0::*)(_Types..., ...) const> { typedef true_type _Bool_type; typedef _Ret result_type; typedef _Arg0 _Class_type; }; template<class _Ret, class _Arg0, class... _Types> struct _Is_memfunptr<_Ret (_Arg0::*)(_Types..., ...) volatile> { typedef true_type _Bool_type; typedef _Ret result_type; typedef _Arg0 _Class_type; }; template<class _Ret, class _Arg0, class... _Types> struct _Is_memfunptr<_Ret (_Arg0::*)(_Types..., ...) const volatile> { typedef true_type _Bool_type; typedef _Ret result_type; typedef _Arg0 _Class_type; }; template<class _Ret, class _Arg0, class... _Types> struct _Is_memfunptr<_Ret (_Arg0::*)(_Types..., ...) &> { typedef true_type _Bool_type; typedef _Ret result_type; typedef _Arg0 _Class_type; }; template<class _Ret, class _Arg0, class... _Types> struct _Is_memfunptr<_Ret (_Arg0::*)(_Types..., ...) const &> { typedef true_type _Bool_type; typedef _Ret result_type; typedef _Arg0 _Class_type; }; template<class _Ret, class _Arg0, class... _Types> struct _Is_memfunptr<_Ret (_Arg0::*)(_Types..., ...) volatile &> { typedef true_type _Bool_type; typedef _Ret result_type; typedef _Arg0 _Class_type; }; template<class _Ret, class _Arg0, class... _Types> struct _Is_memfunptr<_Ret (_Arg0::*)(_Types..., ...) const volatile &> { typedef true_type _Bool_type; typedef _Ret result_type; typedef _Arg0 _Class_type; }; template<class _Ret, class _Arg0, class... _Types> struct _Is_memfunptr<_Ret (_Arg0::*)(_Types..., ...) &&> { typedef true_type _Bool_type; typedef _Ret result_type; typedef _Arg0 _Class_type; }; template<class _Ret, class _Arg0, class... _Types> struct _Is_memfunptr<_Ret (_Arg0::*)(_Types..., ...) const &&> { typedef true_type _Bool_type; typedef _Ret result_type; typedef _Arg0 _Class_type; }; template<class _Ret, class _Arg0, class... _Types> struct _Is_memfunptr<_Ret (_Arg0::*)(_Types..., ...) volatile &&> { typedef true_type _Bool_type; typedef _Ret result_type; typedef _Arg0 _Class_type; }; template<class _Ret, class _Arg0, class... _Types> struct _Is_memfunptr<_Ret (_Arg0::*)(_Types..., ...) const volatile &&> { typedef true_type _Bool_type; typedef _Ret result_type; typedef _Arg0 _Class_type; };


	
template<class _Ty>
	struct is_void
		: false_type
	{	
	};








template<> struct is_void< void> : true_type { }; template<> struct is_void<const void> : true_type { }; template<> struct is_void<volatile void> : true_type { }; template<> struct is_void<const volatile void> : true_type { };


	
	
template<class _Ty>
	struct add_const
	{	
	typedef const _Ty type;
	};

	
template<class _Ty>
	struct add_volatile
	{	
	typedef volatile _Ty type;
	};

	
template<class _Ty>
	struct add_cv
	{	
	typedef const volatile _Ty type;
	};

	
template<class _Ty,
	bool = _Is_function<_Ty>::_Weird || is_void<_Ty>::value>
	struct _Add_reference
	{	
	typedef _Ty _Lvalue;
	typedef _Ty _Rvalue;
	};

template<class _Ty>
	struct _Add_reference<_Ty, false>
	{	
	typedef _Ty& _Lvalue;
	typedef _Ty&& _Rvalue;
	};

	
template<class _Ty>
	struct add_lvalue_reference
	{	
	typedef typename _Add_reference<_Ty>::_Lvalue type;
	};

	
template<class _Ty>
	struct add_rvalue_reference
	{	
	typedef typename _Add_reference<_Ty>::_Rvalue type;
	};

	
template<class _Ty>
	typename add_rvalue_reference<_Ty>::type
		declval() noexcept;

	
template<class _Ty>
	struct remove_extent
	{	
	typedef _Ty type;
	};

template<class _Ty, size_t _Ix>
	struct remove_extent<_Ty[_Ix]>
	{	
	typedef _Ty type;
	};

template<class _Ty>
	struct remove_extent<_Ty[]>
	{	
	typedef _Ty type;
	};

	
template<class _Ty>
	struct remove_all_extents
	{	
	typedef _Ty type;
	};

template<class _Ty, size_t _Ix>
	struct remove_all_extents<_Ty[_Ix]>
	{	
	typedef typename remove_all_extents<_Ty>::type type;
	};

template<class _Ty>
	struct remove_all_extents<_Ty[]>
	{	
	typedef typename remove_all_extents<_Ty>::type type;
	};

	
template<class _Ty>
	struct remove_pointer
	{	
	typedef _Ty type;
	};








template<class _Ty> struct remove_pointer<_Ty *> { typedef _Ty type; }; template<class _Ty> struct remove_pointer<_Ty *const> { typedef _Ty type; }; template<class _Ty> struct remove_pointer<_Ty *volatile> { typedef _Ty type; }; template<class _Ty> struct remove_pointer<_Ty *const volatile> { typedef _Ty type; };


	
template<class _Ty,
	bool = _Is_function<_Ty>::_Weird>
	struct _Add_pointer
	{	
	typedef _Ty type;
	};

template<class _Ty>
	struct _Add_pointer<_Ty, false>
	{	
	typedef typename remove_reference<_Ty>::type *type;
	};

template<class _Ty>
	struct add_pointer
	{	
	typedef typename _Add_pointer<_Ty>::type type;
	};

	
	
template<class _Ty>
	struct is_array
		: false_type
	{	
	};

template<class _Ty, size_t _Nx>
	struct is_array<_Ty[_Nx]>
		: true_type
	{	
	};

template<class _Ty>
	struct is_array<_Ty[]>
		: true_type
	{	
	};

	
template<class _Ty>
	struct is_lvalue_reference
		: false_type
	{	
	};

template<class _Ty>
	struct is_lvalue_reference<_Ty&>
		: true_type
	{	
	};

	
template<class _Ty>
	struct is_rvalue_reference
		: false_type
	{	
	};

template<class _Ty>
	struct is_rvalue_reference<_Ty&&>
		: true_type
	{	
	};

	
template<class _Ty>
	struct is_reference
		: _Cat_base<is_lvalue_reference<_Ty>::value
		|| is_rvalue_reference<_Ty>::value>
	{	
	};


	
template<class _Ty,
	bool _Pmf = _Is_memfunptr<_Ty>::_Bool_type::value>
	struct _Is_member_object_pointer
		: false_type
	{	
	};

template<class _Ty1,
	class _Ty2>
	struct _Is_member_object_pointer<_Ty1 _Ty2::*, false>
		: true_type
	{	
	typedef _Ty2 _Class_type;
	};

template<class _Ty>
	struct is_member_object_pointer
		: _Is_member_object_pointer<typename remove_cv<_Ty>::type>::type
	{	
	};

	
template<class _Ty>
	struct is_member_function_pointer
		: _Is_memfunptr<typename remove_cv<_Ty>::type>::_Bool_type
	{	
	};

	
template<class _Ty>
	struct _Is_pointer
		: false_type
	{	
	};

template<class _Ty>
	struct _Is_pointer<_Ty *>
		: _Cat_base<!is_member_object_pointer<_Ty *>::value
		&& !is_member_function_pointer<_Ty *>::value>
	{	
	};

template<class _Ty>
	struct is_pointer
		: _Is_pointer<typename remove_cv<_Ty>::type>
	{	
	};

	

template<class _Ty>
	struct is_null_pointer
		: _Cat_base<is_same<typename remove_cv<_Ty>::type, nullptr_t>::value>
	{	
	};

	
template<class _Ty>
	struct is_union
		: _Cat_base<__is_union(_Ty)>
	{	
	};

	
template<class _Ty>
	struct is_class
		: _Cat_base<__is_class(_Ty)>
	{	
	};

	
template<class _Ty>
	struct is_fundamental
		: _Cat_base<is_arithmetic<_Ty>::value
		|| is_void<_Ty>::value
		|| is_null_pointer<_Ty>::value>
	{	
	};

	
template<class _Ty>
	struct is_object
		: _Cat_base<!is_function<_Ty>::value
		&& !is_reference<_Ty>::value
		&& !is_void<_Ty>::value>
	{	
	};

	

template<class _From,
	class _To>
	struct is_convertible
		: _Cat_base<__is_convertible_to(_From, _To)>
	{	
	};

	
template<class _Ty>
	struct is_enum
		: _Cat_base<__is_enum(_Ty)>
	{	
	};

	
template<class _Ty>
	struct is_compound
		: _Cat_base<!is_fundamental<_Ty>::value>
	{	
	};

	
template<class _Ty>
	struct is_member_pointer
		: _Cat_base<is_member_object_pointer<_Ty>::value
		|| is_member_function_pointer<_Ty>::value>
	{	
	};

	
template<class _Ty>
	struct is_scalar
		: _Cat_base<is_arithmetic<_Ty>::value
		|| is_enum<_Ty>::value
		|| is_pointer<_Ty>::value
		|| is_member_pointer<_Ty>::value
		|| is_null_pointer<_Ty>::value>
	{	
	};

	
template<class _Ty>
	struct is_const
		: false_type
	{	
	};

template<class _Ty>
	struct is_const<const _Ty>
		: true_type
	{	
	};

	
template<class _Ty>
	struct is_volatile
		: false_type
	{	
	};

template<class _Ty>
	struct is_volatile<volatile _Ty>
		: true_type
	{	
	};

	
template<class _Ty>
	struct is_pod
		: _Cat_base<__is_pod(_Ty)>
	{	
	};

	
template<class _Ty>
	struct is_empty
		: _Cat_base<__is_empty(_Ty)>
	{	
	};

	
template<class _Ty>
	struct is_polymorphic
		: _Cat_base<__is_polymorphic(_Ty)>
	{	
	};

	
template<class _Ty>
	struct is_abstract
		: _Cat_base<__is_abstract(_Ty)>
	{	
	};

	
template<class _Ty>
	struct is_final
		: _Cat_base<__is_final(_Ty)>
	{	
	};

	
template<class _Ty>
	struct is_standard_layout
		: _Cat_base<__is_standard_layout(_Ty)>
	{	
	};

	
template<class _Ty>
	struct is_literal_type
		: _Cat_base<__is_literal_type(_Ty)>
	{	
	};

	
template<class _Ty>
	struct is_trivial
		: _Cat_base<__is_trivial(_Ty)>
	{	
	};

	
template<class _Ty>
	struct is_trivially_copyable
		: _Cat_base<__is_trivially_copyable(_Ty)>
	{	
	};

	
template<class _Ty>
	struct has_virtual_destructor
		: _Cat_base<__has_virtual_destructor(_Ty)>
	{	
	};

		
	

template<class _Ty,
	class... _Args>
	struct is_constructible
		: _Cat_base<__is_constructible(_Ty, _Args...)>
	{	
	};

	
template<class _Ty>
	struct is_copy_constructible
		: is_constructible<
			_Ty,
			typename add_lvalue_reference<
				typename add_const<_Ty>::type
			>::type
		>::type
	{	
	};

	
template<class _Ty>
	struct is_default_constructible
		: is_constructible<_Ty>::type
	{	
	};

	
template<class _Ty>
	struct is_move_constructible
		: is_constructible<
			_Ty,
			typename add_rvalue_reference<_Ty>::type
		>::type
	{	
	};

	
template<class _To,
	class _From>
	struct is_assignable
		: integral_constant<bool, __is_assignable(_To, _From)>
	{	
	};

	
template<class _Ty>
	struct is_copy_assignable
		: is_assignable<
			typename add_lvalue_reference<_Ty>::type,
			typename add_lvalue_reference<
				typename add_const<_Ty>::type
			>::type
		>::type
	{	
	};

	
template<class _Ty>
	struct is_move_assignable
		: is_assignable<
			typename add_lvalue_reference<_Ty>::type,
			typename add_rvalue_reference<_Ty>::type
		>::type
	{	
	};

	
template<class _Ty>
	struct is_destructible
		: _Cat_base<__is_destructible(_Ty)>
	{	
	};

		
	

template<class _Ty,
	class... _Args>
	struct is_trivially_constructible
		: _Cat_base<__is_trivially_constructible(_Ty, _Args...)>
	{	
	};

	
template<class _Ty>
	struct is_trivially_copy_constructible
		: is_trivially_constructible<
			_Ty,
			typename add_lvalue_reference<
				typename add_const<_Ty>::type
			>::type
		>::type
	{	
	};

	
template<class _Ty>
	struct is_trivially_default_constructible
		: is_trivially_constructible<_Ty>::type
	{	
	};

	
template<class _Ty>
	struct is_trivially_move_constructible
		: is_trivially_constructible<
			_Ty,
			typename add_rvalue_reference<_Ty>::type
		>::type
	{	
	};

	
template<class _To,
	class _From>
	struct is_trivially_assignable
		: _Cat_base<__is_trivially_assignable(_To, _From)>
	{	
	};

	
template<class _Ty>
	struct is_trivially_copy_assignable
		: is_trivially_assignable<
			typename add_lvalue_reference<_Ty>::type,
			typename add_lvalue_reference<
				typename add_const<_Ty>::type
			>::type
		>::type
	{	
	};

	
template<class _Ty>
	struct is_trivially_move_assignable
		: is_trivially_assignable<
			typename add_lvalue_reference<_Ty>::type,
			typename add_rvalue_reference<_Ty>::type
		>::type
	{	
	};

	
template<class _Ty>
	struct is_trivially_destructible
		: _Cat_base<__has_trivial_destructor(_Ty)>
	{	
	};

		
	

template<class _Ty,
	class... _Args>
	struct is_nothrow_constructible
		: _Cat_base<__is_nothrow_constructible(_Ty, _Args...)>
	{	
	};

	
template<class _Ty>
	struct is_nothrow_copy_constructible
		: is_nothrow_constructible<
			_Ty,
			typename add_lvalue_reference<
				typename add_const<_Ty>::type
			>::type
		>::type
	{	
	};

	
template<class _Ty>
	struct is_nothrow_default_constructible
		: is_nothrow_constructible<_Ty>::type
	{	
	};

	
template<class _Ty>
	struct is_nothrow_move_constructible
		: is_nothrow_constructible<
			_Ty,
			typename add_rvalue_reference<_Ty>::type
		>::type
	{	
	};

	
template<class _To,
	class _From>
	struct is_nothrow_assignable
		: _Cat_base<__is_nothrow_assignable(_To, _From)>
	{	
	};

	
template<class _Ty>
	struct is_nothrow_copy_assignable
		: is_nothrow_assignable<
			typename add_lvalue_reference<_Ty>::type,
			typename add_lvalue_reference<
				typename add_const<_Ty>::type
			>::type
		>::type
	{	
	};

	
template<class _Ty>
	struct is_nothrow_move_assignable
		: is_nothrow_assignable<
			typename add_lvalue_reference<_Ty>::type,
			typename add_rvalue_reference<_Ty>::type
		>::type
	{	
	};

	
template<class _Ty>
	struct is_nothrow_destructible
		: _Cat_base<__is_nothrow_destructible(_Ty)>
	{	
	};

	
template<class _Ty,
	bool = is_integral<_Ty>::value>
	struct _Sign_base
	{	
	typedef typename remove_cv<_Ty>::type _Uty;
	typedef _Cat_base<_Uty(-1) < _Uty(0)> _Signed;
	typedef _Cat_base<_Uty(0) < _Uty(-1)> _Unsigned;
	};

template<class _Ty>
	struct _Sign_base<_Ty, false>
	{	
		
	typedef is_floating_point<_Ty> _Signed;
	typedef false_type _Unsigned;
	};

template<class _Ty>
	struct is_signed
		: _Sign_base<_Ty>::_Signed
	{	
	};

	
template<class _Ty>
	struct is_unsigned
		: _Sign_base<_Ty>::_Unsigned
	{	
	};

	
template<class _Ty>
	struct _Change_sign
	{	
	static_assert(
		((is_integral<_Ty>::value || is_enum<_Ty>::value)
			&& !is_same<_Ty, bool>::value),
		"make_signed<T>/make_unsigned<T> require that T shall be a (possibly "
		"cv-qualified) integral type or enumeration but not a bool type.");

	typedef
		typename _If<is_same<_Ty, signed char>::value
			|| is_same<_Ty, unsigned char     >::value, signed char,
		typename _If<is_same<_Ty, short       >::value
			|| is_same<_Ty, unsigned short    >::value, short,
		typename _If<is_same<_Ty, int         >::value
			|| is_same<_Ty, unsigned int      >::value, int,
		typename _If<is_same<_Ty, long        >::value
			|| is_same<_Ty, unsigned long     >::value, long,
		typename _If<is_same<_Ty, long long   >::value
			|| is_same<_Ty, unsigned long long>::value, long long,
		typename _If<sizeof (_Ty) == sizeof (signed char), signed char,
		typename _If<sizeof (_Ty) == sizeof (short      ), short,
		typename _If<sizeof (_Ty) == sizeof (int        ), int,
		typename _If<sizeof (_Ty) == sizeof (long       ), long,
			long long
		>::type>::type>::type>::type>::type>::type>::type>::type>::type
			_Signed;

	typedef
		typename _If<is_same<_Signed, signed char>::value, unsigned char,
		typename _If<is_same<_Signed, short      >::value, unsigned short,
		typename _If<is_same<_Signed, int        >::value, unsigned int,
		typename _If<is_same<_Signed, long       >::value, unsigned long,
			unsigned long long
		>::type>::type>::type>::type
			_Unsigned;
	};

template<class _Ty>
	struct _Change_sign<const _Ty>
	{	
	typedef const typename _Change_sign<_Ty>::_Signed _Signed;
	typedef const typename _Change_sign<_Ty>::_Unsigned _Unsigned;
	};

template<class _Ty>
	struct _Change_sign<volatile _Ty>
	{	
	typedef volatile typename _Change_sign<_Ty>::_Signed _Signed;
	typedef volatile typename _Change_sign<_Ty>::_Unsigned _Unsigned;
	};

template<class _Ty>
	struct _Change_sign<const volatile _Ty>
	{	
	typedef const volatile typename _Change_sign<_Ty>::_Signed _Signed;
	typedef const volatile typename _Change_sign<_Ty>::_Unsigned _Unsigned;
	};

	
template<class _Ty>
	struct make_signed
	{	
	typedef typename _Change_sign<_Ty>::_Signed type;
	};

	
template<class _Ty>
	struct make_unsigned
	{	
	typedef typename _Change_sign<_Ty>::_Unsigned type;
	};

	

template<class _Ty>
	struct alignment_of
		: integral_constant<size_t, alignof(_Ty)>
	{	
	};

	




template<class _Ty,
	size_t _Len>
	union _Align_type
	{	
	_Ty _Val;
	char _Pad[_Len];
	};

template<size_t _Len,
	size_t _Align,
	class _Ty,
	bool _Ok>
	struct _Aligned;

template<size_t _Len,
	size_t _Align,
	class _Ty>
	struct _Aligned<_Len, _Align, _Ty, true>
	{	
	typedef _Align_type<_Ty, _Len> type;
	};

template<size_t _Len,
	size_t _Align>
	struct _Aligned<_Len, _Align, double, false>
	{	
	typedef _Align_type<max_align_t, _Len> type;
	};

template<size_t _Len,
	size_t _Align>
	struct _Aligned<_Len, _Align, int, false>
	{	
	typedef typename _Aligned<_Len, _Align, double, _Align <= alignment_of<double>::value>::type type;
	};

template<size_t _Len,
	size_t _Align>
	struct _Aligned<_Len, _Align, short, false>
	{	
	typedef typename _Aligned<_Len, _Align, int, _Align <= alignment_of<int>::value>::type type;
	};

template<size_t _Len,
	size_t _Align>
	struct _Aligned<_Len, _Align, char, false>
	{	
	typedef typename _Aligned<_Len, _Align, short, _Align <= alignment_of<short>::value>::type type;
	};

template<size_t _Len,
	size_t _Align = alignment_of<max_align_t>::value>
	struct aligned_storage
	{	
	typedef typename _Aligned<_Len, _Align, char, _Align <= alignment_of<char>::value>::type type;
	};




	
template<size_t... _Vals>
	struct _Maximum;

template<>
	struct _Maximum<>
	{	
	static constexpr size_t value = 0;
	};

template<size_t _Val>
	struct _Maximum<_Val>
	{	
	static constexpr size_t value = _Val;
	};

template<size_t _First,
	size_t _Second,
	size_t... _Rest>
	struct _Maximum<_First, _Second, _Rest...>
		: _Maximum<(_First < _Second ? _Second : _First), _Rest...>
	{	
	};

template<size_t _Len,
	class... _Types>
	struct aligned_union
	{	
	static constexpr size_t _Max_len = _Maximum<
		_Len, sizeof(_Types)...>::value;	
	static constexpr size_t alignment_value = _Maximum<
		alignment_of<_Types>::value...>::value;
	typedef typename aligned_storage<_Max_len, alignment_value>::type type;
	};

	
template<class _Ty>
	struct underlying_type
	{	
	typedef __underlying_type(_Ty) type;
	};

	
template<class _Ty>
	struct rank
		: integral_constant<size_t, 0>
	{	
	};

template<class _Ty, size_t _Ix>
	struct rank<_Ty[_Ix]>
		: integral_constant<size_t, rank<_Ty>::value + 1>
	{	
	};

template<class _Ty>
	struct rank<_Ty[]>
		: integral_constant<size_t, rank<_Ty>::value + 1>
	{	
	};

	
template<class _Ty, unsigned int _Nx>
	struct _Extent
		: integral_constant<size_t, 0>
	{	
	};

template<class _Ty, size_t _Ix>
	struct _Extent<_Ty[_Ix], 0>
		: integral_constant<size_t, _Ix>
	{	
	};

template<class _Ty, unsigned int _Nx, size_t _Ix>
	struct _Extent<_Ty[_Ix], _Nx>
		: _Extent<_Ty, _Nx - 1>
	{	
	};

template<class _Ty, unsigned int _Nx>
	struct _Extent<_Ty[], _Nx>
		: _Extent<_Ty, _Nx - 1>
	{	
	};

template<class _Ty, unsigned int _Nx = 0>
	struct extent
		: _Extent<_Ty, _Nx>
	{	
	};

	
template<class _Base,
	class _Der>
	struct is_base_of
		: _Cat_base<__is_base_of(_Base, _Der)>
	{	
	};

	
template<class _Ty>
	struct decay
	{	
	typedef typename remove_reference<_Ty>::type _Ty1;

	typedef typename _If<is_array<_Ty1>::value,
		typename remove_extent<_Ty1>::type *,
		typename _If<is_function<_Ty1>::value,
			typename add_pointer<_Ty1>::type,
			typename remove_cv<_Ty1>::type>::type>::type type;
	};

	
template<class...>
	struct _Conjunction;

template<bool,
	class _Lhs,
	class... _Traits>
	struct _Choose_conjunction
	{	
	typedef _Lhs type;
	};

template<class _Lhs,
	class... _Traits>
	struct _Choose_conjunction<true, _Lhs, _Traits...>
	{	
	typedef typename _Conjunction<_Traits...>::type type;
	};

template<>
	struct _Conjunction<>
	{	
	typedef true_type type;
	};

template<class _Trait>
	struct _Conjunction<_Trait>
	{	
	typedef _Trait type;
	};

template<class _Lhs,
	class... _Traits>
	struct _Conjunction<_Lhs, _Traits...>
	{	
	typedef typename _Choose_conjunction<_Lhs::value, _Lhs, _Traits...>::type type;
	};

template<class... _Traits>
	struct conjunction
		: _Conjunction<_Traits...>::type
	{	
		
		
	};

	
template<class...>
	struct _Disjunction;

template<bool,
	class _Lhs,
	class... _Traits>
	struct _Choose_disjunction
	{	
	typedef _Lhs type;
	};

template<class _Lhs,
	class... _Traits>
	struct _Choose_disjunction<false, _Lhs, _Traits...>
	{	
	typedef typename _Disjunction<_Traits...>::type type;
	};

template<>
	struct _Disjunction<>
	{	
	typedef false_type type;
	};

template<class _Trait>
	struct _Disjunction<_Trait>
	{	
	typedef _Trait type;
	};

template<class _Lhs,
	class... _Traits>
	struct _Disjunction<_Lhs, _Traits...>
	{	
	typedef typename _Choose_disjunction<_Lhs::value, _Lhs, _Traits...>::type type;
	};

template<class... _Traits>
	struct disjunction
		: _Disjunction<_Traits...>::type
	{	
		
		
	};

	
template<class _Trait>
	struct negation
		: bool_constant<!_Trait::value>
	{	
	};


namespace tr1 {	
using ::std:: add_const;
using ::std:: add_cv;
using ::std:: add_pointer;
using ::std:: add_volatile;
using ::std:: aligned_storage;
using ::std:: alignment_of;
using ::std:: conditional;
using ::std:: decay;
using ::std:: enable_if;
using ::std:: extent;
using ::std:: false_type;
using ::std:: has_virtual_destructor;
using ::std:: integral_constant;
using ::std:: is_abstract;
using ::std:: is_arithmetic;
using ::std:: is_array;
using ::std:: is_base_of;
using ::std:: is_class;
using ::std:: is_compound;
using ::std:: is_const;
using ::std:: is_convertible;
using ::std:: is_empty;
using ::std:: is_enum;
using ::std:: is_floating_point;
using ::std:: is_function;
using ::std:: is_fundamental;
using ::std:: is_integral;
using ::std:: is_member_function_pointer;
using ::std:: is_member_object_pointer;
using ::std:: is_member_pointer;
using ::std:: is_object;
using ::std:: is_pod;
using ::std:: is_pointer;
using ::std:: is_polymorphic;
using ::std:: is_reference;
using ::std:: is_same;
using ::std:: is_scalar;
using ::std:: is_signed;
using ::std:: is_union;
using ::std:: is_unsigned;
using ::std:: is_void;
using ::std:: is_volatile;
using ::std:: make_signed;
using ::std:: make_unsigned;
using ::std:: rank;
using ::std:: remove_all_extents;
using ::std:: remove_const;
using ::std:: remove_cv;
using ::std:: remove_extent;
using ::std:: remove_pointer;
using ::std:: remove_reference;
using ::std:: remove_volatile;
using ::std:: true_type;
	}	


		
template<class... _Ty>
	struct common_type;

template<class _Ty>
	struct common_type<_Ty>
	{	
	typedef typename decay<_Ty>::type type;
	};

template<class _Ty0,
	class _Ty1>
	struct common_type<_Ty0, _Ty1>
	{	
	typedef typename decay<
		decltype(_Always_false<_Ty0>::value
			? ::std:: declval<_Ty0>()
			: ::std:: declval<_Ty1>())
	>::type type;
	};

template<class _Ty0,
	class _Ty1,
	class... _Ty>
	struct common_type<_Ty0, _Ty1, _Ty...>
	{	
	typedef typename common_type<
		typename common_type<_Ty0, _Ty1>::type, _Ty...
	>::type type;
	};

	
template<class _Ty,
	_Ty... _Vals>
	struct integer_sequence
	{	
	static_assert(is_integral<_Ty>::value,
		"integer_sequence<T, I...> requires T to be an integral type.");

	typedef integer_sequence<_Ty, _Vals...> type;
	typedef _Ty value_type;

	static constexpr size_t size() noexcept
		{	
		return (sizeof...(_Vals));
		}
	};

	
 




































template<class _Ty,
	_Ty _Size>
	using make_integer_sequence = __make_integer_seq<integer_sequence, _Ty, _Size>;
 

template<size_t... _Vals>
	using index_sequence = integer_sequence<size_t, _Vals...>;

template<size_t _Size>
	using make_index_sequence = make_integer_sequence<size_t, _Size>;

template<class... _Types>
	using index_sequence_for = make_index_sequence<sizeof...(_Types)>;


	
template<class _Ty>
	struct identity
	{	
	typedef _Ty type;

	const _Ty& operator()(const _Ty& _Left) const
		{	
		return (_Left);
		}
	};


	
template<class _Ty> inline
	constexpr _Ty&& forward(
		typename remove_reference<_Ty>::type& _Arg) noexcept
	{	
	return (static_cast<_Ty&&>(_Arg));
	}

template<class _Ty> inline
	constexpr _Ty&& forward(
		typename remove_reference<_Ty>::type&& _Arg) noexcept
	{	
	static_assert(!is_lvalue_reference<_Ty>::value, "bad forward call");
	return (static_cast<_Ty&&>(_Arg));
	}

		
template<class _Ty> inline
	constexpr typename remove_reference<_Ty>::type&&
		move(_Ty&& _Arg) noexcept
	{	
	return (static_cast<typename remove_reference<_Ty>::type&&>(_Arg));
	}

		
template<class _Ty> inline
	constexpr typename _If<!is_nothrow_move_constructible<_Ty>::value
		&& is_copy_constructible<_Ty>::value,
			const _Ty&, _Ty&&>::type
	move_if_noexcept(_Ty& _Arg) noexcept
	{	
	return (::std:: move(_Arg));
	}

	
template<class...>
	struct _Param_tester
	{	
	typedef void type;
	};

	
template<class... _Types>	
	using void_t = typename _Param_tester<_Types...>::type;

	
struct _Invoker_pmf_object
	{	
	template<class _Decayed,
		class _Ty1,
		class... _Types2>
		static auto _Call(_Decayed _Pmf, _Ty1&& _Arg1, _Types2&&... _Args2)
		-> decltype((::std:: forward<_Ty1>(_Arg1).*_Pmf)(
			::std:: forward<_Types2>(_Args2)...))
		{	
		return ((::std:: forward<_Ty1>(_Arg1).*_Pmf)(
			::std:: forward<_Types2>(_Args2)...));
		}
	};

struct _Invoker_pmf_pointer
	{	
	template<class _Decayed,
		class _Ty1,
		class... _Types2>
		static auto _Call(_Decayed _Pmf, _Ty1&& _Arg1, _Types2&&... _Args2)
		-> decltype(((*::std:: forward<_Ty1>(_Arg1)).*_Pmf)(
			::std:: forward<_Types2>(_Args2)...))
		{	
		return (((*::std:: forward<_Ty1>(_Arg1)).*_Pmf)(
			::std:: forward<_Types2>(_Args2)...));
		}
	};

struct _Invoker_pmd_object
	{	
	template<class _Decayed,
		class _Ty1>
		static auto _Call(_Decayed _Pmd, _Ty1&& _Arg1)
		-> decltype(::std:: forward<_Ty1>(_Arg1).*_Pmd)
		{	
		return (::std:: forward<_Ty1>(_Arg1).*_Pmd);
		}
	};

struct _Invoker_pmd_pointer
	{	
	template<class _Decayed,
		class _Ty1>
		static auto _Call(_Decayed _Pmd, _Ty1&& _Arg1)
		-> decltype((*::std:: forward<_Ty1>(_Arg1)).*_Pmd)
		{	
		return ((*::std:: forward<_Ty1>(_Arg1)).*_Pmd);
		}
	};

struct _Invoker_functor
	{	
	template<class _Callable,
		class... _Types>
		static auto _Call(_Callable&& _Obj, _Types&&... _Args)
		-> decltype(::std:: forward<_Callable>(_Obj)(
			::std:: forward<_Types>(_Args)...))
		{	
		return (::std:: forward<_Callable>(_Obj)(
			::std:: forward<_Types>(_Args)...));
		}
	};

template<class _Callable,
	class _Ty1,
	class _Decayed = typename decay<_Callable>::type,
	bool _Is_pmf = is_member_function_pointer<_Decayed>::value,
	bool _Is_pmd = is_member_object_pointer<_Decayed>::value>
	struct _Invoker1;

template<class _Callable,
	class _Ty1,
	class _Decayed>
	struct _Invoker1<_Callable, _Ty1, _Decayed, true, false>
		: _If<is_base_of<
			typename _Is_memfunptr<_Decayed>::_Class_type,
			typename decay<_Ty1>::type>::value,
		_Invoker_pmf_object,
		_Invoker_pmf_pointer>::type
	{	
	};

template<class _Callable,
	class _Ty1,
	class _Decayed>
	struct _Invoker1<_Callable, _Ty1, _Decayed, false, true>
		: _If<is_base_of<
			typename _Is_member_object_pointer<_Decayed>::_Class_type,
			typename decay<_Ty1>::type>::value,
		_Invoker_pmd_object,
		_Invoker_pmd_pointer>::type
	{	
	};

template<class _Callable,
	class _Ty1,
	class _Decayed>
	struct _Invoker1<_Callable, _Ty1, _Decayed, false, false>
		: _Invoker_functor
	{	
	};

template<class _Callable,
	class... _Types>
	struct _Invoker;

template<class _Callable>
	struct _Invoker<_Callable>
		: _Invoker_functor
	{	
	};

template<class _Callable,
	class _Ty1,
	class... _Types2>
	struct _Invoker<_Callable, _Ty1, _Types2...>
		: _Invoker1<_Callable, _Ty1>
	{	
	};

template<class _Callable,
	class... _Types> inline
	auto invoke(_Callable&& _Obj, _Types&&... _Args)
	-> decltype(_Invoker<_Callable, _Types...>::_Call(
		::std:: forward<_Callable>(_Obj), ::std:: forward<_Types>(_Args)...))
	{	
	return (_Invoker<_Callable, _Types...>::_Call(
		::std:: forward<_Callable>(_Obj), ::std:: forward<_Types>(_Args)...));
	}

template<class _Rx,
	bool = is_void<_Rx>::value>
	struct _Forced
	{	
	};

struct _Unforced
	{	
	};

template<class _Cv_void,
	class... _Valtys> inline
	void _Invoke_ret(_Forced<_Cv_void, true>, _Valtys&&... _Vals)
	{	
	::std:: invoke(::std:: forward<_Valtys>(_Vals)...);
	}

template<class _Rx,
	class... _Valtys> inline
	_Rx _Invoke_ret(_Forced<_Rx, false>, _Valtys&&... _Vals)
	{	
	return (::std:: invoke(::std:: forward<_Valtys>(_Vals)...));
	}

template<class... _Valtys> inline
	auto _Invoke_ret(_Forced<_Unforced, false>, _Valtys&&... _Vals)
	-> decltype(::std:: invoke(::std:: forward<_Valtys>(_Vals)...))
	{	
	return (::std:: invoke(::std:: forward<_Valtys>(_Vals)...));
	}

	
struct _Unique_tag_result_of
	{	
	};

template<class _Void,
	class... _Types>
	struct _Result_of
	{	
	};

template<class... _Types>
	struct _Result_of<
		void_t<
			_Unique_tag_result_of,	
			decltype(::std:: invoke(::std:: declval<_Types>()...))>,
		_Types...>
	{	
	typedef decltype(::std:: invoke(::std:: declval<_Types>()...)) type;
	};

template<class _Fty>
	struct result_of
	{	
	static_assert(_Always_false<_Fty>::value,
		"result_of<CallableType> is invalid; use "
		"result_of<CallableType(zero or more argument types)> instead.");
	};









template<class _Fty, class... _Args> struct result_of<_Fty __cdecl (_Args...)> : _Result_of<void, _Fty, _Args...> { };    template<class _Fty, class... _Args> struct result_of<_Fty __vectorcall (_Args...)> : _Result_of<void, _Fty, _Args...> { };


	
template<class _Ty,
	class = void>
	struct _Weak_result_type
	{	
	};

template<class _Ty>
	struct _Weak_result_type<_Ty, void_t<
		typename _Ty::result_type> >
	{	
	typedef typename _Ty::result_type result_type;
	};

template<class _Ty,
	class = void>
	struct _Weak_argument_type
		: _Weak_result_type<_Ty>
	{	
	};

template<class _Ty>
	struct _Weak_argument_type<_Ty, void_t<
		typename _Ty::argument_type> >
		: _Weak_result_type<_Ty>
	{	
	typedef typename _Ty::argument_type argument_type;
	};

template<class _Ty,
	class = void>
	struct _Weak_binary_args
		: _Weak_argument_type<_Ty>
	{	
	};

template<class _Ty>
	struct _Weak_binary_args<_Ty, void_t<
		typename _Ty::first_argument_type,
		typename _Ty::second_argument_type> >
		: _Weak_argument_type<_Ty>
	{	
	typedef typename _Ty::first_argument_type first_argument_type;
	typedef typename _Ty::second_argument_type second_argument_type;
	};

template<class _Ty>
	struct _Weak_types
	{	
	typedef _Is_function<typename remove_pointer<_Ty>::type> _Is_f_or_pf;
	typedef _Is_memfunptr<typename remove_cv<_Ty>::type> _Is_pmf;
	typedef typename _If<_Is_f_or_pf::_Bool_type::value, _Is_f_or_pf,
		typename _If<_Is_pmf::_Bool_type::value, _Is_pmf,
		_Weak_binary_args<_Ty> >::type>::type type;
	};

	
template<class _Ty>
	class reference_wrapper
		: public _Weak_types<_Ty>::type
	{	
public:
	static_assert(is_object<_Ty>::value || is_function<_Ty>::value,
		"reference_wrapper<T> requires T to be an object type "
		"or a function type.");

	typedef _Ty type;

	reference_wrapper(_Ty& _Val) noexcept
		: _Ptr(::std:: addressof(_Val))
		{	
		}

	operator _Ty&() const noexcept
		{	
		return (*_Ptr);
		}

	_Ty& get() const noexcept
		{	
		return (*_Ptr);
		}

	template<class... _Types>
		auto operator()(_Types&&... _Args) const
		-> decltype(::std:: invoke(get(), ::std:: forward<_Types>(_Args)...))
		{	
		return (::std:: invoke(get(), ::std:: forward<_Types>(_Args)...));
		}

	reference_wrapper(_Ty&&) = delete;

private:
	_Ty *_Ptr;
	};

	
template<class _Ty> inline
	reference_wrapper<_Ty>
		ref(_Ty& _Val) noexcept
	{	
	return (reference_wrapper<_Ty>(_Val));
	}

template<class _Ty>
	void ref(const _Ty&&) = delete;

template<class _Ty> inline
	reference_wrapper<_Ty>
		ref(reference_wrapper<_Ty> _Val) noexcept
	{	
	return (::std:: ref(_Val.get()));
	}

template<class _Ty> inline
	reference_wrapper<const _Ty>
		cref(const _Ty& _Val) noexcept
	{	
	return (reference_wrapper<const _Ty>(_Val));
	}

template<class _Ty>
	void cref(const _Ty&&) = delete;

template<class _Ty> inline
	reference_wrapper<const _Ty>
		cref(reference_wrapper<_Ty> _Val) noexcept
	{	
	return (::std:: cref(_Val.get()));
	}

	
template<class _Ty>
	struct _Unrefwrap_helper
	{	
	typedef _Ty type;
	static constexpr bool _Is_refwrap = false;
	};

template<class _Ty>
	struct _Unrefwrap_helper<reference_wrapper<_Ty> >
	{	
	typedef _Ty& type;
	static constexpr bool _Is_refwrap = true;
	};

template<class _Ty>
	struct _Unrefwrap
	{	
	typedef typename decay<_Ty>::type _Ty1;
	typedef typename _Unrefwrap_helper<_Ty1>::type type;
	static constexpr bool _Is_refwrap = _Unrefwrap_helper<_Ty1>::_Is_refwrap;
	};


namespace tr1 {	
using ::std:: cref;
using ::std:: ref;
using ::std:: reference_wrapper;
using ::std:: result_of;
	}	


		
template<class _Ty>
	struct _Is_swappable;

		
template<class _Ty>
	struct _Is_nothrow_swappable;

		





template<class _Ty,
	class = void> inline

	void swap(_Ty&, _Ty&)
		noexcept(is_nothrow_move_constructible<_Ty>::value && is_nothrow_move_assignable<_Ty>::value);

template<class _Ty,
	size_t _Size,
	class = typename enable_if<_Is_swappable<_Ty>::value>::type> inline
	void swap(_Ty (&)[_Size], _Ty (&)[_Size])
		noexcept(_Is_nothrow_swappable<_Ty>::value);

		
template<class _Ty1,
	class _Ty2,
	class = void>
	struct _Swappable_with_helper
		: false_type
	{	
	};

struct _Swappable_with_helper_unique_type {}; 
template<class _Ty1,
	class _Ty2>
	struct _Swappable_with_helper<_Ty1, _Ty2, void_t<
		_Swappable_with_helper_unique_type,
		decltype(swap(::std:: declval<_Ty1>(), ::std:: declval<_Ty2>()))>>
		: true_type
	{	
	};

		
template<class _Ty1,
	class _Ty2>
	struct _Is_swappable_with
		: conjunction<
			_Swappable_with_helper<_Ty1, _Ty2>,
			_Swappable_with_helper<_Ty2, _Ty1>>::type
	{	
		
	};

		
template<class _Ty>
	struct _Is_swappable
		: _Is_swappable_with<
			typename add_lvalue_reference<_Ty>::type,
			typename add_lvalue_reference<_Ty>::type>::type
	{	
	};

		
template<class _Ty1,
	class _Ty2>
	struct _Swap_cannot_throw
	{	
		
		

	static constexpr bool value = 
		noexcept(swap(::std:: declval<_Ty1>(), ::std:: declval<_Ty2>()))
		&& noexcept(swap(::std:: declval<_Ty2>(), ::std:: declval<_Ty1>()));



	using type = bool_constant<value>;
	};

		
template<class _Ty1,
	class _Ty2>
	struct _Is_nothrow_swappable_with
		: conjunction<
			_Is_swappable_with<_Ty1, _Ty2>,
			_Swap_cannot_throw<_Ty1, _Ty2>>::type
	{	
		
	};

		
template<class _Ty>
	struct _Is_nothrow_swappable
		: _Is_nothrow_swappable_with<
			typename add_lvalue_reference<_Ty>::type,
			typename add_lvalue_reference<_Ty>::type>::type
	{	
	};



































		
template<class _Ty>
	using remove_const_t = typename remove_const<_Ty>::type;

template<class _Ty>
	using remove_volatile_t = typename remove_volatile<_Ty>::type;

template<class _Ty>
	using remove_cv_t = typename remove_cv<_Ty>::type;

template<class _Ty>
	using add_const_t = typename add_const<_Ty>::type;

template<class _Ty>
	using add_volatile_t = typename add_volatile<_Ty>::type;

template<class _Ty>
	using add_cv_t = typename add_cv<_Ty>::type;

template<class _Ty>
	using remove_reference_t = typename remove_reference<_Ty>::type;

template<class _Ty>
	using add_lvalue_reference_t = typename add_lvalue_reference<_Ty>::type;

template<class _Ty>
	using add_rvalue_reference_t = typename add_rvalue_reference<_Ty>::type;

template<class _Ty>
	using make_signed_t = typename make_signed<_Ty>::type;

template<class _Ty>
	using make_unsigned_t = typename make_unsigned<_Ty>::type;

template<class _Ty>
	using remove_extent_t = typename remove_extent<_Ty>::type;

template<class _Ty>
	using remove_all_extents_t = typename remove_all_extents<_Ty>::type;

template<class _Ty>
	using remove_pointer_t = typename remove_pointer<_Ty>::type;

template<class _Ty>
	using add_pointer_t = typename add_pointer<_Ty>::type;

template<size_t _Len,
	size_t _Align = alignment_of<max_align_t>::value>
	using aligned_storage_t = typename aligned_storage<_Len, _Align>::type;

template<size_t _Len,
	class... _Types>
	using aligned_union_t = typename aligned_union<_Len, _Types...>::type;

template<class _Ty>
	using decay_t = typename decay<_Ty>::type;

template<bool _Test,
	class _Ty = void>
	using enable_if_t = typename enable_if<_Test, _Ty>::type;

template<bool _Test,
	class _Ty1,
	class _Ty2>
	using conditional_t = typename conditional<_Test, _Ty1, _Ty2>::type;

template<class... _Ty>
	using common_type_t = typename common_type<_Ty...>::type;

template<class _Ty>
	using underlying_type_t = typename underlying_type<_Ty>::type;

template<class _Ty>
	using result_of_t = typename result_of<_Ty>::type;

	
 
template<class _Ty>
	constexpr bool is_void_v = is_void<_Ty>::value;
template<class _Ty>
	constexpr bool is_null_pointer_v = is_null_pointer<_Ty>::value;
template<class _Ty>
	constexpr bool is_array_v = is_array<_Ty>::value;
template<class _Ty>
	constexpr bool is_pointer_v = is_pointer<_Ty>::value;
template<class _Ty>
	constexpr bool is_lvalue_reference_v = is_lvalue_reference<_Ty>::value;
template<class _Ty>
	constexpr bool is_rvalue_reference_v = is_rvalue_reference<_Ty>::value;
template<class _Ty>
	constexpr bool is_member_object_pointer_v = is_member_object_pointer<_Ty>::value;
template<class _Ty>
	constexpr bool is_member_function_pointer_v = is_member_function_pointer<_Ty>::value;
template<class _Ty>
	constexpr bool is_enum_v = is_enum<_Ty>::value;
template<class _Ty>
	constexpr bool is_union_v = is_union<_Ty>::value;
template<class _Ty>
	constexpr bool is_class_v = is_class<_Ty>::value;
template<class _Ty>
	constexpr bool is_reference_v = is_reference<_Ty>::value;
template<class _Ty>
	constexpr bool is_fundamental_v = is_fundamental<_Ty>::value;
template<class _Ty>
	constexpr bool is_object_v = is_object<_Ty>::value;
template<class _Ty>
	constexpr bool is_scalar_v = is_scalar<_Ty>::value;
template<class _Ty>
	constexpr bool is_compound_v = is_compound<_Ty>::value;
template<class _Ty>
	constexpr bool is_member_pointer_v = is_member_pointer<_Ty>::value;
template<class _Ty>
	constexpr bool is_const_v = is_const<_Ty>::value;
template<class _Ty>
	constexpr bool is_volatile_v = is_volatile<_Ty>::value;
template<class _Ty>
	constexpr bool is_trivial_v = is_trivial<_Ty>::value;
template<class _Ty>
	constexpr bool is_trivially_copyable_v = is_trivially_copyable<_Ty>::value;
template<class _Ty>
	constexpr bool is_standard_layout_v = is_standard_layout<_Ty>::value;
template<class _Ty>
	constexpr bool is_pod_v = is_pod<_Ty>::value;
template<class _Ty>
	constexpr bool is_literal_type_v = is_literal_type<_Ty>::value;
template<class _Ty>
	constexpr bool is_empty_v = is_empty<_Ty>::value;
template<class _Ty>
	constexpr bool is_polymorphic_v = is_polymorphic<_Ty>::value;
template<class _Ty>
	constexpr bool is_abstract_v = is_abstract<_Ty>::value;
template<class _Ty>
	constexpr bool is_final_v = is_final<_Ty>::value;
template<class _Ty>
	constexpr bool is_signed_v = is_signed<_Ty>::value;
template<class _Ty>
	constexpr bool is_unsigned_v = is_unsigned<_Ty>::value;
template<class _Ty,
	class... _Args>
	constexpr bool is_constructible_v = is_constructible<_Ty, _Args...>::value;
template<class _Ty>
	constexpr bool is_default_constructible_v = is_default_constructible<_Ty>::value;
template<class _Ty>
	constexpr bool is_copy_constructible_v = is_copy_constructible<_Ty>::value;
template<class _Ty>
	constexpr bool is_move_constructible_v = is_move_constructible<_Ty>::value;
template<class _Ty,
	class _Uty>
	constexpr bool is_assignable_v = is_assignable<_Ty, _Uty>::value;
template<class _Ty>
	constexpr bool is_copy_assignable_v = is_copy_assignable<_Ty>::value;
template<class _Ty>
	constexpr bool is_move_assignable_v = is_move_assignable<_Ty>::value;







template<class _Ty>
	constexpr bool is_destructible_v = is_destructible<_Ty>::value;
template<class _Ty,
	class... _Args>
	constexpr bool is_trivially_constructible_v = is_trivially_constructible<_Ty, _Args...>::value;
template<class _Ty>
	constexpr bool is_trivially_default_constructible_v = is_trivially_default_constructible<_Ty>::value;
template<class _Ty>
	constexpr bool is_trivially_copy_constructible_v = is_trivially_copy_constructible<_Ty>::value;
template<class _Ty>
	constexpr bool is_trivially_move_constructible_v = is_trivially_move_constructible<_Ty>::value;
template<class _Ty,
	class _Uty>
	constexpr bool is_trivially_assignable_v = is_trivially_assignable<_Ty, _Uty>::value;
template<class _Ty>
	constexpr bool is_trivially_copy_assignable_v = is_trivially_copy_assignable<_Ty>::value;
template<class _Ty>
	constexpr bool is_trivially_move_assignable_v = is_trivially_move_assignable<_Ty>::value;
template<class _Ty>
	constexpr bool is_trivially_destructible_v = is_trivially_destructible<_Ty>::value;
template<class _Ty,
	class... _Args>
	constexpr bool is_nothrow_constructible_v = is_nothrow_constructible<_Ty, _Args...>::value;
template<class _Ty>
	constexpr bool is_nothrow_default_constructible_v = is_nothrow_default_constructible<_Ty>::value;
template<class _Ty>
	constexpr bool is_nothrow_copy_constructible_v = is_nothrow_copy_constructible<_Ty>::value;
template<class _Ty>
	constexpr bool is_nothrow_move_constructible_v = is_nothrow_move_constructible<_Ty>::value;
template<class _Ty,
	class _Uty>
	constexpr bool is_nothrow_assignable_v = is_nothrow_assignable<_Ty, _Uty>::value;
template<class _Ty>
	constexpr bool is_nothrow_copy_assignable_v = is_nothrow_copy_assignable<_Ty>::value;
template<class _Ty>
	constexpr bool is_nothrow_move_assignable_v = is_nothrow_move_assignable<_Ty>::value;







template<class _Ty>
	constexpr bool is_nothrow_destructible_v = is_nothrow_destructible<_Ty>::value;
template<class _Ty>
	constexpr bool has_virtual_destructor_v = has_virtual_destructor<_Ty>::value;
template<class _Ty>
	constexpr size_t alignment_of_v = alignment_of<_Ty>::value;
template<class _Ty>
	constexpr size_t rank_v = rank<_Ty>::value;
template<class _Ty,
	unsigned int _Ix = 0>
	constexpr size_t extent_v = extent<_Ty, _Ix>::value;
template<class _Base,
	class _Derived>
	constexpr bool is_base_of_v = is_base_of<_Base, _Derived>::value;
template<class _From,
	class _To>
	constexpr bool is_convertible_v = is_convertible<_From, _To>::value;
template<class... _Traits>
	constexpr bool conjunction_v = conjunction<_Traits...>::value;
template<class... _Traits>
	constexpr bool disjunction_v = disjunction<_Traits...>::value;
template<class _Trait>
	constexpr bool negation_v = negation<_Trait>::value;
 

}

 
 #pragma warning(pop)
 #pragma pack(pop)










 #pragma pack(push,8)
 #pragma warning(push,3)
 
 

namespace std {

  


  



  




  


  

__declspec(dllimport) bool __cdecl uncaught_exception() noexcept;
__declspec(dllimport) int __cdecl uncaught_exceptions() noexcept;

}

 









#pragma once





__pragma(pack(push, 8)) extern "C" {





    


















typedef struct _heapinfo
{
    int* _pentry;
    size_t _size;
    int _useflag;
} _HEAPINFO;








   
void* __cdecl _alloca(  size_t _Size);





    __declspec(dllimport) intptr_t __cdecl _get_heap_handle(void);

     
    __declspec(dllimport) int __cdecl _heapmin(void);

    
        __declspec(dllimport) int __cdecl _heapwalk(  _HEAPINFO* _EntryInfo);
    

    
          __declspec(dllimport) int __cdecl _heapchk(void);
        __declspec(dllimport) int __cdecl _resetstkoflw(void);
    
     
    
    
    

    
        
    



    typedef char __static_assert_t[(sizeof(unsigned int) <= 16) != 0];


    #pragma warning(push)
    #pragma warning(disable:6540)

    __inline void* _MarkAllocaS(   void* _Ptr, unsigned int _Marker)
    {
        if (_Ptr)
        {
            *((unsigned int*)_Ptr) = _Marker;
            _Ptr = (char*)_Ptr + 16;
        }
        return _Ptr;
    }

    __inline size_t _MallocaComputeSize(size_t _Size)
    {
        size_t _MarkedSize = _Size + 16;
        return _MarkedSize > _Size ? _MarkedSize : 0;
    }

    #pragma warning(pop)


















    
    















    

    #pragma warning(push)
    #pragma warning(disable: 6014)
    __inline void __cdecl _freea(    void* _Memory)
    {
        unsigned int _Marker;
        if (_Memory)
        {
            _Memory = (char*)_Memory - 16;
            _Marker = *(unsigned int*)_Memory;
            if (_Marker == 0xDDDD)
            {
                free(_Memory);
            }
            





        }
    }
    #pragma warning(pop)






    




} __pragma(pack(pop))









#pragma once









#pragma once










#pragma once





__pragma(pack(push, 8)) extern "C" {



typedef void (__cdecl* terminate_handler )();
typedef void (__cdecl* terminate_function)();








    __declspec(dllimport) __declspec(noreturn) void __cdecl abort();
    __declspec(dllimport) __declspec(noreturn) void __cdecl terminate() throw();

    

        __declspec(dllimport) terminate_handler __cdecl set_terminate(
              terminate_handler _NewTerminateHandler
            ) throw();

        __declspec(dllimport) terminate_handler __cdecl _get_terminate();

    



} __pragma(pack(pop))






__pragma(pack(push, 8)) extern "C" {



typedef void (__cdecl* unexpected_handler )();
typedef void (__cdecl* unexpected_function)();






struct _EXCEPTION_POINTERS;


    
    __declspec(dllimport) __declspec(noreturn) void __cdecl unexpected() throw(...);

    

        __declspec(dllimport) unexpected_handler __cdecl set_unexpected(
              unexpected_handler _NewUnexpectedHandler
            ) throw();

        __declspec(dllimport) unexpected_handler __cdecl _get_unexpected();

        typedef void (__cdecl* _se_translator_function)(unsigned int, struct _EXCEPTION_POINTERS*);

        __declspec(dllimport) _se_translator_function __cdecl _set_se_translator(
              _se_translator_function _NewSETranslator
            );

    

    class type_info;

    __declspec(dllimport) int __cdecl _is_exception_typeof(
          type_info const&     _Type,
          _EXCEPTION_POINTERS* _ExceptionPtr
        );

    __declspec(dllimport) bool __cdecl __uncaught_exception();
    __declspec(dllimport) int  __cdecl __uncaught_exceptions();



} __pragma(pack(pop))








#pragma pack(push, 8)


__pragma(pack(push, 8)) extern "C" {

struct __std_exception_data
{
    char const* _What;
    bool        _DoFree;
};

__declspec(dllimport) void __cdecl __std_exception_copy(
       __std_exception_data const* _From,
      __std_exception_data*       _To
    );

__declspec(dllimport) void __cdecl __std_exception_destroy(
      __std_exception_data* _Data
    );

} __pragma(pack(pop))



namespace std {

class exception
{
public:

    exception() throw()
        : _Data()
    {
    }

    explicit exception(char const* const _Message) throw()
        : _Data()
    {
        __std_exception_data _InitData = { _Message, true };
        __std_exception_copy(&_InitData, &_Data);
    }

    exception(char const* const _Message, int) throw()
        : _Data()
    {
        _Data._What = _Message;
    }

    exception(exception const& _Other) throw()
        : _Data()
    {
        __std_exception_copy(&_Other._Data, &_Data);
    }

    exception& operator=(exception const& _Other) throw()
    {
        if (this == &_Other)
        {
            return *this;
        }

        __std_exception_destroy(&_Data);
        __std_exception_copy(&_Other._Data, &_Data);
        return *this;
    }

    virtual ~exception() throw()
    {
        __std_exception_destroy(&_Data);
    }

    virtual char const* what() const
    {
        return _Data._What ? _Data._What : "Unknown exception";
    }

private:

    __std_exception_data _Data;
};

class bad_exception
    : public exception
{
public:

    bad_exception() throw()
        : exception("bad exception", 1)
    {
    }
};

class bad_alloc
    : public exception
{
public:

    bad_alloc() throw()
        : exception("bad allocation", 1)
    {
    }

private:

    friend class bad_array_new_length;

    bad_alloc(char const* const _Message) throw()
        : exception(_Message, 1)
    {
    }
};

class bad_array_new_length
    : public bad_alloc
{
public:

    bad_array_new_length() throw()
        : bad_alloc("bad array new length")
    {
    }
};

} 


#pragma pack(pop)







namespace std {

using ::set_terminate; using ::terminate_handler; using ::terminate; using ::set_unexpected; using ::unexpected_handler; using ::unexpected;

typedef void (__cdecl *_Prhand)(const exception&);


inline terminate_handler __cdecl get_terminate() noexcept
	{	
	return (_get_terminate());
	}

inline unexpected_handler __cdecl get_unexpected() noexcept
	{	
	return (_get_unexpected());
	}


}

 
















































































































































































__declspec(dllimport) void __cdecl __ExceptionPtrCreate(  void*);
__declspec(dllimport) void __cdecl __ExceptionPtrDestroy(  void*);
__declspec(dllimport) void __cdecl __ExceptionPtrCopy(  void*,   const void*);
__declspec(dllimport) void __cdecl __ExceptionPtrAssign(  void*,   const void*);
__declspec(dllimport) bool __cdecl __ExceptionPtrCompare(  const void*,   const void*);
__declspec(dllimport) bool __cdecl __ExceptionPtrToBool(  const void*);
__declspec(dllimport) void __cdecl __ExceptionPtrSwap(  void*,   void*);
__declspec(dllimport) void __cdecl __ExceptionPtrCurrentException(  void*);
[[noreturn]] __declspec(dllimport) void __cdecl __ExceptionPtrRethrow(  const void*);
__declspec(dllimport) void __cdecl __ExceptionPtrCopyException(
	  void*,   const void*,   const void*);

namespace std {

class exception_ptr
	{
public:
	exception_ptr() throw ()
		{
		__ExceptionPtrCreate(this);
		}

	exception_ptr(nullptr_t) throw ()
		{
		__ExceptionPtrCreate(this);
		}

	~exception_ptr() throw ()
		{
		__ExceptionPtrDestroy(this);
		}

	exception_ptr(const exception_ptr& _Rhs) throw ()
		{
		__ExceptionPtrCopy(this, &_Rhs);
		}

	exception_ptr& operator=(const exception_ptr& _Rhs) throw ()
		{
		__ExceptionPtrAssign(this, &_Rhs);
		return *this;
		}

	exception_ptr& operator=(nullptr_t) throw ()
		{
		exception_ptr _Ptr;
		__ExceptionPtrAssign(this, &_Ptr);
		return *this;
		}

	explicit operator bool() const throw ()
		{
		return __ExceptionPtrToBool(this);
		}

	[[noreturn]] void _RethrowException() const
		{
		__ExceptionPtrRethrow(this);
		}

	static exception_ptr _Current_exception() throw ()
		{
		exception_ptr _Retval;
		__ExceptionPtrCurrentException(&_Retval);
		return _Retval;
		}

	static exception_ptr _Copy_exception(  void* _Except,   const void* _Ptr)
		{
		exception_ptr _Retval = 0;
		if (!_Ptr)
			{
			
			return _Retval;
			}
		__ExceptionPtrCopyException(&_Retval, _Except, _Ptr);
		return _Retval;
		}

private:
	void* _Data1;
	void* _Data2;
	};

inline void swap(exception_ptr& _Lhs, exception_ptr& _Rhs) throw ()
	{
	__ExceptionPtrSwap(&_Lhs, &_Rhs);
	}

inline bool operator==(const exception_ptr& _Lhs, const exception_ptr& _Rhs) throw ()
	{
	return __ExceptionPtrCompare(&_Lhs, &_Rhs);
	}

inline bool operator==(nullptr_t, const exception_ptr& _Rhs) throw ()
	{
	return !_Rhs;
	}

inline bool operator==(const exception_ptr& _Lhs, nullptr_t) throw ()
	{
	return !_Lhs;
	}

inline bool operator!=(const exception_ptr& _Lhs, const exception_ptr& _Rhs) throw ()
	{
	return !(_Lhs == _Rhs);
	}

inline bool operator!=(nullptr_t _Lhs, const exception_ptr& _Rhs) throw ()
	{
	return !(_Lhs == _Rhs);
	}

inline bool operator!=(const exception_ptr& _Lhs, nullptr_t _Rhs) throw ()
	{
	return !(_Lhs == _Rhs);
	}

inline exception_ptr current_exception() noexcept
	{
	return exception_ptr::_Current_exception();
	}

[[noreturn]] inline void rethrow_exception(  exception_ptr _Ptr)
	{
	_Ptr._RethrowException();
	}

template<class _Ex> void *__GetExceptionInfo(_Ex);

template<class _Ex> exception_ptr make_exception_ptr(_Ex _Except) noexcept
	{
	return exception_ptr::_Copy_exception(::std:: addressof(_Except), __GetExceptionInfo(_Except));
	}

	
class nested_exception
	{	
public:
	nested_exception() noexcept
		: _Exc(::std:: current_exception())
		{	
		}

	nested_exception(const nested_exception&) noexcept = default;
	nested_exception& operator=(const nested_exception&) noexcept = default;
	virtual ~nested_exception() noexcept = default;

	[[noreturn]] void rethrow_nested() const
		{	
		if (_Exc)
			::std:: rethrow_exception(_Exc);
		else
			::std:: terminate();
		}

	::std:: exception_ptr nested_ptr() const noexcept
		{	
		return (_Exc);
		}

private:
	::std:: exception_ptr _Exc;
	};

	
template<class _Ty,
	class _Uty>
	struct _With_nested
		: _Uty, nested_exception
	{	
	explicit _With_nested(_Ty&& _Arg)
		: _Uty(::std:: forward<_Ty>(_Arg)), nested_exception()
		{	
		}
	};

template<class _Ty>
	[[noreturn]] inline void _Throw_with_nested(_Ty&& _Arg, true_type)
	{	
	typedef typename remove_reference<_Ty>::type _Uty;
	typedef _With_nested<_Ty, _Uty> _Glued;

	throw _Glued(::std:: forward<_Ty>(_Arg));
	}

template<class _Ty>
	[[noreturn]] inline void _Throw_with_nested(_Ty&& _Arg, false_type)
	{	
	typedef typename decay<_Ty>::type _Decayed;

	throw _Decayed(::std:: forward<_Ty>(_Arg));
	}

template<class _Ty>
	[[noreturn]] inline void throw_with_nested(_Ty&& _Arg)
	{	
	typedef typename remove_reference<_Ty>::type _Uty;

	integral_constant<bool,
		is_class<_Uty>::value
		&& !is_base_of<nested_exception, _Uty>::value
		&& !is_final<_Uty>::value> _Tag;

	_Throw_with_nested(::std:: forward<_Ty>(_Arg), _Tag);
	}

	
template<class _Ty> inline
	void _Rethrow_if_nested(const _Ty *_Ptr, true_type)
	{	
	const auto _Nested = dynamic_cast<const nested_exception *>(_Ptr);

	if (_Nested)
		_Nested->rethrow_nested();
	}

template<class _Ty> inline
	void _Rethrow_if_nested(const _Ty *, false_type)
	{	
	}

template<class _Ty> inline
	void rethrow_if_nested(const _Ty& _Arg)
	{	
	integral_constant<bool,
		is_polymorphic<_Ty>::value
		&& (!is_base_of<nested_exception, _Ty>::value
			|| is_convertible<_Ty *, nested_exception *>::value)> _Tag;

	_Rethrow_if_nested(::std:: addressof(_Arg), _Tag);
	}
}

 
 #pragma warning(pop)
 #pragma pack(pop)

















#pragma once




































































































































































































































































































































extern "C++" {

#pragma pack(push, 8)

#pragma warning(push)
#pragma warning(disable: 4985) 






    namespace std
    {
        struct nothrow_t { };

        extern nothrow_t const nothrow;
    }


   
__declspec(allocator) void* __cdecl operator new(
    size_t _Size
    );

     
__declspec(allocator) void* __cdecl operator new(
    size_t                _Size,
    std::nothrow_t const&
    ) throw();

   
__declspec(allocator) void* __cdecl operator new[](
    size_t _Size
    );

     
__declspec(allocator) void* __cdecl operator new[](
    size_t                _Size,
    std::nothrow_t const&
    ) throw();

void __cdecl operator delete(
    void* _Block
    ) throw();

void __cdecl operator delete(
    void* _Block,
    std::nothrow_t const&
    ) throw();

void __cdecl operator delete[](
    void* _Block
    ) throw();

void __cdecl operator delete[](
    void* _Block,
    std::nothrow_t const&
    ) throw();

void __cdecl operator delete(
    void*  _Block,
    size_t _Size
    ) throw();

void __cdecl operator delete[](
    void* _Block,
    size_t _Size
    ) throw();


    
       
    inline void* __cdecl operator new(size_t _Size,   void* _Where) throw()
    {
        (void)_Size;
        return _Where;
    }

    inline void __cdecl operator delete(void*, void*) throw()
    {
        return;
    }



    
       
    inline void* __cdecl operator new[](size_t _Size,   void* _Where) throw()
    {
        (void)_Size;
        return _Where;
    }

    inline void __cdecl operator delete[](void*, void*) throw()
    {
    }




#pragma warning(pop)
#pragma pack(pop)

} 



 #pragma pack(push,8)
 #pragma warning(push,3)
 

  



namespace std {

		
 

typedef void (__cdecl * new_handler) ();
 

		
__declspec(dllimport) new_handler __cdecl set_new_handler(  new_handler)
	noexcept;	

__declspec(dllimport) new_handler __cdecl get_new_handler()
	noexcept;	
}

 
 #pragma warning(pop)
 #pragma pack(pop)











#pragma once







#pragma once






#pragma once





#pragma once










 







#pragma once





__pragma(pack(push, 8)) extern "C" {































    
















    
    



typedef __int64 fpos_t;




__declspec(dllimport) errno_t __cdecl _get_stream_buffer_pointers(
           FILE*   _Stream,
      char*** _Base,
      char*** _Pointer,
      int**   _Count
    );









    
    __declspec(dllimport) errno_t __cdecl clearerr_s(
          FILE* _Stream
        );

    
    __declspec(dllimport) errno_t __cdecl fopen_s(
          FILE**      _Stream,
                             char const* _FileName,
                             char const* _Mode
        );
    
    
     
    __declspec(dllimport) size_t __cdecl fread_s(
            void*  _Buffer,
                       size_t _BufferSize,
                                                                        size_t _ElementSize,
                                                                        size_t _ElementCount,
                                                                     FILE*  _Stream
        );
    
    
    __declspec(dllimport) errno_t __cdecl freopen_s(
          FILE**      _Stream,
                             char const* _FileName,
                             char const* _Mode,
                            FILE*       _OldStream
        );

     
    __declspec(dllimport) char* __cdecl gets_s(
          char*   _Buffer,
                           rsize_t _Size
        );

    
    __declspec(dllimport) errno_t __cdecl tmpfile_s(
            FILE** _Stream
        );

     
    
    __declspec(dllimport) errno_t __cdecl tmpnam_s(
          char*   _Buffer,
                           rsize_t _Size
        );



__declspec(dllimport) void __cdecl clearerr(
      FILE* _Stream
    );

 

__declspec(dllimport) int __cdecl fclose(
      FILE* _Stream
    );


__declspec(dllimport) int __cdecl _fcloseall(void);

 
__declspec(dllimport) FILE* __cdecl _fdopen(
        int         _FileHandle,
      char const* _Mode
    );

 
__declspec(dllimport) int __cdecl feof(
      FILE* _Stream
    );

 
__declspec(dllimport) int __cdecl ferror(
      FILE* _Stream
    );


__declspec(dllimport) int __cdecl fflush(
      FILE* _Stream
    );

 

__declspec(dllimport) int __cdecl fgetc(
      FILE* _Stream
    );


__declspec(dllimport) int __cdecl _fgetchar(void);

 

__declspec(dllimport) int __cdecl fgetpos(
      FILE*   _Stream,
        fpos_t* _Position
    );

 

__declspec(dllimport) char* __cdecl fgets(
      char* _Buffer,
                           int   _MaxCount,
                        FILE* _Stream
    );

 
__declspec(dllimport) int __cdecl _fileno(
      FILE* _Stream
    );


__declspec(dllimport) int __cdecl _flushall(void);

  __declspec(deprecated("This function or variable may be unsafe. Consider using " "fopen_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__declspec(dllimport) FILE* __cdecl fopen(
      char const* _FileName,
      char const* _Mode
    );


 

__declspec(dllimport) int __cdecl fputc(
         int   _Character,
      FILE* _Stream
    );


__declspec(dllimport) int __cdecl _fputchar(
      int _Character
    );

 

__declspec(dllimport) int __cdecl fputs(
       char const* _Buffer,
      FILE*       _Stream
    );


__declspec(dllimport) size_t __cdecl fread(
      void*  _Buffer,
                                                  size_t _ElementSize,
                                                  size_t _ElementCount,
                                               FILE*  _Stream
    );

 
  __declspec(deprecated("This function or variable may be unsafe. Consider using " "freopen_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__declspec(dllimport) FILE* __cdecl freopen(
       char const* _FileName,
       char const* _Mode,
      FILE*       _Stream
    );

 
__declspec(dllimport) FILE* __cdecl _fsopen(
      char const* _FileName,
      char const* _Mode,
        int         _ShFlag
    );

 

__declspec(dllimport) int __cdecl fsetpos(
      FILE*         _Stream,
         fpos_t const* _Position
    );

 

__declspec(dllimport) int __cdecl fseek(
      FILE* _Stream,
         long  _Offset,
         int   _Origin
    );

 

__declspec(dllimport) int __cdecl _fseeki64(
      FILE*   _Stream,
         __int64 _Offset,
         int     _Origin
    );

 
 
__declspec(dllimport) long __cdecl ftell(
      FILE* _Stream
    );

 
 
__declspec(dllimport) __int64 __cdecl _ftelli64(
      FILE* _Stream
    );


__declspec(dllimport) size_t __cdecl fwrite(
      void const* _Buffer,
                                                size_t      _ElementSize,
                                                size_t      _ElementCount,
                                             FILE*       _Stream
    );

 
 
__declspec(dllimport) int __cdecl getc(
      FILE* _Stream
    );

 
__declspec(dllimport) int __cdecl getchar(void);

 
__declspec(dllimport) int __cdecl _getmaxstdio(void);

extern "C++" { template <size_t _Size> inline char* __cdecl gets_s(char (&_Buffer)[_Size]) throw() { return gets_s(_Buffer, _Size); } }

 
__declspec(dllimport) int __cdecl _getw(
      FILE* _Stream
    );

__declspec(dllimport) void __cdecl perror(
      char const* _ErrorMessage
    );



     
    
    __declspec(dllimport) int __cdecl _pclose(
          FILE* _Stream
        );

     
    __declspec(dllimport) FILE* __cdecl _popen(
          char const* _Command,
          char const* _Mode
        );



 

__declspec(dllimport) int __cdecl putc(
         int   _Character,
      FILE* _Stream
    );


__declspec(dllimport) int __cdecl putchar(
      int _Character
    );


__declspec(dllimport) int __cdecl puts(
      char const* _Buffer
    );

 

__declspec(dllimport) int __cdecl _putw(
         int   _Word, 
      FILE* _Stream
    );



__declspec(dllimport) int __cdecl remove(
      char const* _FileName
    );

 
__declspec(dllimport) int __cdecl rename(
      char const* _OldFileName,
      char const* _NewFileName
    );

__declspec(dllimport) int __cdecl _unlink(
      char const* _FileName
    );



    __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_unlink" ". See online help for details."))
    __declspec(dllimport) int __cdecl unlink(
          char const* _FileName
        );





__declspec(dllimport) void __cdecl rewind(
      FILE* _Stream
    );


__declspec(dllimport) int __cdecl _rmtmp(void);

__declspec(deprecated("This function or variable may be unsafe. Consider using " "setvbuf" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__declspec(dllimport) void __cdecl setbuf(
                                                  FILE* _Stream,
        char* _Buffer
    );


__declspec(dllimport) int __cdecl _setmaxstdio(
      int _Maximum
    );

 

__declspec(dllimport) int __cdecl setvbuf(
                           FILE*  _Stream,
        char*  _Buffer,
                              int    _Mode,
                              size_t _Size
    );






 
__declspec(dllimport) __declspec(allocator) char* __cdecl _tempnam(
      char const* _DirectoryName,
      char const* _FilePrefix
    );





  __declspec(deprecated("This function or variable may be unsafe. Consider using " "tmpfile_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__declspec(dllimport) FILE* __cdecl tmpfile(void);

extern "C++" { template <size_t _Size> inline errno_t __cdecl tmpnam_s(  char (&_Buffer)[_Size]) throw() { return tmpnam_s(_Buffer, _Size); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "tmpnam_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport)  char* __cdecl tmpnam(  char *_Buffer);

 

__declspec(dllimport) int __cdecl ungetc(
         int   _Character,
      FILE* _Stream
    );








__declspec(dllimport) void __cdecl _lock_file(
      FILE* _Stream
    );

__declspec(dllimport) void __cdecl _unlock_file(
      FILE* _Stream
    );

 

__declspec(dllimport) int __cdecl _fclose_nolock(
      FILE* _Stream
    );

 

__declspec(dllimport) int __cdecl _fflush_nolock(
      FILE* _Stream
    );

 

__declspec(dllimport) int __cdecl _fgetc_nolock(
      FILE* _Stream
    );

 

__declspec(dllimport) int __cdecl _fputc_nolock(
         int   _Character,
      FILE* _Stream
    );


__declspec(dllimport) size_t __cdecl _fread_nolock(
      void*  _Buffer,
                                                  size_t _ElementSize,
                                                  size_t _ElementCount,
                                               FILE*  _Stream
    );


 
__declspec(dllimport) size_t __cdecl _fread_nolock_s(
      void*  _Buffer,
               size_t _BufferSize,
                                                                  size_t _ElementSize,
                                                                  size_t _ElementCount,
                                                               FILE*  _Stream
    );


__declspec(dllimport) int __cdecl _fseek_nolock(
      FILE* _Stream,
         long  _Offset,
         int   _Origin
    );


__declspec(dllimport) int __cdecl _fseeki64_nolock(
      FILE*   _Stream,
         __int64 _Offset,
         int     _Origin
    );

 
__declspec(dllimport) long __cdecl _ftell_nolock(
      FILE* _Stream
    );

 
__declspec(dllimport) __int64 __cdecl _ftelli64_nolock(
      FILE* _Stream
    );


__declspec(dllimport) size_t __cdecl _fwrite_nolock(
      void const* _Buffer,
                                                size_t      _ElementSize,
                                                size_t      _ElementCount,
                                             FILE*       _Stream
    );


__declspec(dllimport) int __cdecl _getc_nolock(
      FILE* _Stream
    );


__declspec(dllimport) int __cdecl _putc_nolock(
         int   _Character,
      FILE* _Stream
    );


__declspec(dllimport) int __cdecl _ungetc_nolock(
         int   _Character,
      FILE* _Stream
    );


























__declspec(dllimport) int* __cdecl __p__commode(void);




    














__declspec(dllimport) int __cdecl __stdio_common_vfprintf(
                                         unsigned __int64 _Options,
                                      FILE*            _Stream,
        char const*      _Format,
                                     _locale_t        _Locale,
                                            va_list          _ArgList
    );

__declspec(dllimport) int __cdecl __stdio_common_vfprintf_s(
                                         unsigned __int64 _Options,
                                      FILE*            _Stream,
        char const*      _Format,
                                     _locale_t        _Locale,
                                            va_list          _ArgList
    );

 
__declspec(dllimport) int __cdecl __stdio_common_vfprintf_p(
                                         unsigned __int64 _Options,
                                      FILE*            _Stream,
        char const*      _Format,
                                     _locale_t        _Locale,
                                            va_list          _ArgList
    );


__inline int __cdecl _vfprintf_l(
       FILE*       const _Stream,
        char const* const _Format,
      _locale_t   const _Locale,
             va_list           _ArgList
    )



{
    return __stdio_common_vfprintf((*__local_stdio_printf_options()), _Stream, _Format, _Locale, _ArgList);
}



__inline int __cdecl vfprintf(
                            FILE*       const _Stream,
        char const* const _Format,
                                  va_list           _ArgList
    )



{
    return _vfprintf_l(_Stream, _Format, 0, _ArgList);
}



__inline int __cdecl _vfprintf_s_l(
       FILE*       const _Stream,
        char const* const _Format,
      _locale_t   const _Locale,
             va_list           _ArgList
    )



{
    return __stdio_common_vfprintf_s((*__local_stdio_printf_options()), _Stream, _Format, _Locale, _ArgList);
}




    
    __inline int __cdecl vfprintf_s(
                                FILE*       const _Stream,
            char const* const _Format,
                                      va_list           _ArgList
        )



    {
        return _vfprintf_s_l(_Stream, _Format, 0, _ArgList);
    }





__inline int __cdecl _vfprintf_p_l(
       FILE*       const _Stream,
        char const* const _Format,
      _locale_t   const _Locale,
             va_list           _ArgList
    )



{
    return __stdio_common_vfprintf_p((*__local_stdio_printf_options()), _Stream, _Format, _Locale, _ArgList);
}



__inline int __cdecl _vfprintf_p(
                            FILE*       const _Stream,
        char const* const _Format,
                                  va_list           _ArgList
    )



{
    return _vfprintf_p_l(_Stream, _Format, 0, _ArgList);
}



__inline int __cdecl _vprintf_l(
        char const* const _Format,
                                     _locale_t   const _Locale,
                                            va_list           _ArgList
    )



{
    return _vfprintf_l((__acrt_iob_func(1)), _Format, _Locale, _ArgList);
}



__inline int __cdecl vprintf(
        char const* const _Format,
                                  va_list           _ArgList
    )



{
    return _vfprintf_l((__acrt_iob_func(1)), _Format, 0, _ArgList);
}



__inline int __cdecl _vprintf_s_l(
        char const* const _Format,
                                     _locale_t   const _Locale,
                                            va_list           _ArgList
    )



{
    return _vfprintf_s_l((__acrt_iob_func(1)), _Format, _Locale, _ArgList);
}




    
    __inline int __cdecl vprintf_s(
            char const* const _Format,
                                      va_list           _ArgList
        )



    {
        return _vfprintf_s_l((__acrt_iob_func(1)), _Format, 0, _ArgList);
    }





__inline int __cdecl _vprintf_p_l(
        char const* const _Format,
                                     _locale_t   const _Locale,
                                            va_list           _ArgList
    )



{
    return _vfprintf_p_l((__acrt_iob_func(1)), _Format, _Locale, _ArgList);
}



__inline int __cdecl _vprintf_p(
        char const* const _Format,
                                  va_list           _ArgList
    )



{
    return _vfprintf_p_l((__acrt_iob_func(1)), _Format, 0, _ArgList);
}



__inline int __cdecl _fprintf_l(
                                      FILE*       const _Stream,
        char const* const _Format,
                                     _locale_t   const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vfprintf_l(_Stream, _Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}



__inline int __cdecl fprintf(
                            FILE*       const _Stream,
        char const* const _Format,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
    _Result = _vfprintf_l(_Stream, _Format, 0, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}


__declspec(dllimport) int __cdecl _set_printf_count_output(
      int _Value
    );

__declspec(dllimport) int __cdecl _get_printf_count_output(void);


__inline int __cdecl _fprintf_s_l(
                                      FILE*       const _Stream,
        char const* const _Format,
                                     _locale_t   const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vfprintf_s_l(_Stream, _Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}




    
    __inline int __cdecl fprintf_s(
                                FILE*       const _Stream,
            char const* const _Format,
        ...)



    {
        int _Result;
        va_list _ArgList;
        ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
        _Result = _vfprintf_s_l(_Stream, _Format, 0, _ArgList);
        ((void)(_ArgList = (va_list)0));
        return _Result;
    }





__inline int __cdecl _fprintf_p_l(
                                      FILE*       const _Stream,
        char const* const _Format,
                                     _locale_t   const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vfprintf_p_l(_Stream, _Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}



__inline int __cdecl _fprintf_p(
                            FILE*       const _Stream,
        char const* const _Format,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
    _Result = _vfprintf_p_l(_Stream, _Format, 0, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}



__inline int __cdecl _printf_l(
        char const* const _Format,
                                     _locale_t   const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vfprintf_l((__acrt_iob_func(1)), _Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}



__inline int __cdecl printf(
        char const* const _Format,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
    _Result = _vfprintf_l((__acrt_iob_func(1)), _Format, 0, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}



__inline int __cdecl _printf_s_l(
        char const* const _Format,
                                     _locale_t   const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vfprintf_s_l((__acrt_iob_func(1)), _Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}




    
    __inline int __cdecl printf_s(
            char const* const _Format,
        ...)



    {
        int _Result;
        va_list _ArgList;
        ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
        _Result = _vfprintf_s_l((__acrt_iob_func(1)), _Format, 0, _ArgList);
        ((void)(_ArgList = (va_list)0));
        return _Result;
    }





__inline int __cdecl _printf_p_l(
        char const* const _Format,
                                     _locale_t   const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vfprintf_p_l((__acrt_iob_func(1)), _Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}



__inline int __cdecl _printf_p(
        char const* const _Format,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
    _Result = _vfprintf_p_l((__acrt_iob_func(1)), _Format, 0, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}








__declspec(dllimport) int __cdecl __stdio_common_vfscanf(
                                        unsigned __int64 _Options,
                                     FILE*            _Stream,
        char const*      _Format,
                                    _locale_t        _Locale,
                                           va_list          _Arglist
    );


__inline int __cdecl _vfscanf_l(
                            FILE*       const _Stream,
        char const* const _Format,
                           _locale_t   const _Locale,
                                  va_list           _ArgList
    )



{
    return __stdio_common_vfscanf(
        (*__local_stdio_scanf_options ()),
        _Stream, _Format, _Locale, _ArgList);
}



__inline int __cdecl vfscanf(
                            FILE*       const _Stream,
        char const* const _Format,
                                  va_list           _ArgList
    )



{
    return _vfscanf_l(_Stream, _Format, 0, _ArgList);
}



__inline int __cdecl _vfscanf_s_l(
                            FILE*       const _Stream,
        char const* const _Format,
                           _locale_t   const _Locale,
                                  va_list           _ArgList
    )



{
    return __stdio_common_vfscanf(
        (*__local_stdio_scanf_options ()) | (1ULL << 0),
        _Stream, _Format, _Locale, _ArgList);
}





    
    __inline int __cdecl vfscanf_s(
                                FILE*       const _Stream,
            char const* const _Format,
                                      va_list           _ArgList
        )



    {
        return _vfscanf_s_l(_Stream, _Format, 0, _ArgList);
    }





__inline int __cdecl _vscanf_l(
        char const* const _Format,
                           _locale_t   const _Locale,
                                  va_list           _ArgList
    )



{
    return _vfscanf_l((__acrt_iob_func(0)), _Format, _Locale, _ArgList);
}



__inline int __cdecl vscanf(
        char const* const _Format,
                                  va_list           _ArgList
    )



{
    return _vfscanf_l((__acrt_iob_func(0)), _Format, 0, _ArgList);
}



__inline int __cdecl _vscanf_s_l(
        char const* const _Format,
                           _locale_t   const _Locale,
                                  va_list           _ArgList
    )



{
    return _vfscanf_s_l((__acrt_iob_func(0)), _Format, _Locale, _ArgList);
}




    
    __inline int __cdecl vscanf_s(
            char const* const _Format,
                                      va_list           _ArgList
        )



    {
        return _vfscanf_s_l((__acrt_iob_func(0)), _Format, 0, _ArgList);
    }




 __declspec(deprecated("This function or variable may be unsafe. Consider using " "_fscanf_s_l" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__inline int __cdecl _fscanf_l(
                                     FILE*       const _Stream,
        char const* const _Format,
                                    _locale_t   const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vfscanf_l(_Stream, _Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}


  __declspec(deprecated("This function or variable may be unsafe. Consider using " "fscanf_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__inline int __cdecl fscanf(
                           FILE*       const _Stream,
        char const* const _Format,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
    _Result = _vfscanf_l(_Stream, _Format, 0, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}



__inline int __cdecl _fscanf_s_l(
                                       FILE*       const _Stream,
        char const* const _Format,
                                      _locale_t   const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vfscanf_s_l(_Stream, _Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}




    
    __inline int __cdecl fscanf_s(
                                 FILE*       const _Stream,
            char const* const _Format,
        ...)



    {
        int _Result;
        va_list _ArgList;
        ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
        _Result = _vfscanf_s_l(_Stream, _Format, 0, _ArgList);
        ((void)(_ArgList = (va_list)0));
        return _Result;
    }




 __declspec(deprecated("This function or variable may be unsafe. Consider using " "_scanf_s_l" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__inline int __cdecl _scanf_l(
        char const* const _Format,
                                    _locale_t   const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vfscanf_l((__acrt_iob_func(0)), _Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}


  __declspec(deprecated("This function or variable may be unsafe. Consider using " "scanf_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__inline int __cdecl scanf(
        char const* const _Format,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
    _Result = _vfscanf_l((__acrt_iob_func(0)), _Format, 0, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}



__inline int __cdecl _scanf_s_l(
        char const* const _Format,
                                      _locale_t   const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vfscanf_s_l((__acrt_iob_func(0)), _Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}




    
    __inline int __cdecl scanf_s(
            char const* const _Format,
        ...)



    {
        int _Result;
        va_list _ArgList;
        ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
        _Result = _vfscanf_s_l((__acrt_iob_func(0)), _Format, 0, _ArgList);
        ((void)(_ArgList = (va_list)0));
        return _Result;
    }











 
__declspec(dllimport) int __cdecl __stdio_common_vsprintf(
                                         unsigned __int64 _Options,
                 char*            _Buffer,
                                         size_t           _BufferCount,
        char const*      _Format,
                                     _locale_t        _Locale,
                                            va_list          _ArgList
    );

 
__declspec(dllimport) int __cdecl __stdio_common_vsprintf_s(
                                         unsigned __int64 _Options,
                 char*            _Buffer,
                                         size_t           _BufferCount,
        char const*      _Format,
                                     _locale_t        _Locale,
                                            va_list          _ArgList
    );

 
__declspec(dllimport) int __cdecl __stdio_common_vsnprintf_s(
                                         unsigned __int64 _Options,
                 char*            _Buffer,
                                         size_t           _BufferCount,
                                         size_t           _MaxCount,
        char const*      _Format,
                                     _locale_t        _Locale,
                                            va_list          _ArgList
    );

 
__declspec(dllimport) int __cdecl __stdio_common_vsprintf_p(
                                         unsigned __int64 _Options,
                 char*            _Buffer,
                                         size_t           _BufferCount,
        char const*      _Format,
                                     _locale_t        _Locale,
                                            va_list          _ArgList
    );

 
 __declspec(deprecated("This function or variable may be unsafe. Consider using " "_vsnprintf_s_l" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__inline int __cdecl _vsnprintf_l(
                   char*       const _Buffer,
                                         size_t      const _BufferCount,
        char const* const _Format,
                                     _locale_t   const _Locale,
                                            va_list           _ArgList
    )



{
    int const _Result = __stdio_common_vsprintf(
        (*__local_stdio_printf_options()) | (1ULL << 0),
        _Buffer, _BufferCount, _Format, _Locale, _ArgList);

    return _Result < 0 ? -1 : _Result;
}


 

__inline int __cdecl _vsnprintf(
        char*       const _Buffer,
                                          size_t      const _BufferCount,
                   char const* const _Format,
                                             va_list           _ArgList
    )



{
    #pragma warning(push)
    #pragma warning(disable: 4996) 
    return _vsnprintf_l(_Buffer, _BufferCount, _Format, 0, _ArgList);
    #pragma warning(pop)
}













 

__inline int __cdecl vsnprintf(
         char*       const _Buffer,
                               size_t      const _BufferCount,
        char const* const _Format,
                                  va_list           _ArgList
    )



{
    int const _Result = __stdio_common_vsprintf(
        (*__local_stdio_printf_options()) | (1ULL << 1),
        _Buffer, _BufferCount, _Format, 0, _ArgList);

    return _Result < 0 ? -1 : _Result;
}


 
 __declspec(deprecated("This function or variable may be unsafe. Consider using " "_vsprintf_s_l" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__inline int __cdecl _vsprintf_l(
        char*       const _Buffer,
                      char const* const _Format,
                    _locale_t   const _Locale,
                           va_list           _ArgList
    )



{
    #pragma warning(push)
    #pragma warning(disable: 4996) 
    return _vsnprintf_l(_Buffer, (size_t)-1, _Format, _Locale, _ArgList);
    #pragma warning(pop)
}


 
 __declspec(deprecated("This function or variable may be unsafe. Consider using " "vsprintf_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__inline int __cdecl vsprintf(
               char*       const _Buffer,
        char const* const _Format,
                                  va_list           _ArgList
    )



{
    #pragma warning(push)
    #pragma warning(disable: 4996) 
    return _vsnprintf_l(_Buffer, (size_t)-1, _Format, 0, _ArgList);
    #pragma warning(pop)
}


 

__inline int __cdecl _vsprintf_s_l(
                 char*       const _Buffer,
                                         size_t      const _BufferCount,
        char const* const _Format,
                                     _locale_t   const _Locale,
                                            va_list           _ArgList
    )



{
    int const _Result = __stdio_common_vsprintf_s(
        (*__local_stdio_printf_options()),
        _Buffer, _BufferCount, _Format, _Locale, _ArgList);

    return _Result < 0 ? -1 : _Result;
}




     
    
    __inline int __cdecl vsprintf_s(
           char*       const _Buffer,
                                   size_t      const _BufferCount,
            char const* const _Format,
                                      va_list           _ArgList
        )



    {
        return _vsprintf_s_l(_Buffer, _BufferCount, _Format, 0, _ArgList);
    }

    
    extern "C++" { template <size_t _Size> inline   int __cdecl vsprintf_s(  char (&_Buffer)[_Size],     char const* _Format, va_list _ArgList) throw() { return vsprintf_s(_Buffer, _Size, _Format, _ArgList); } }



 

__inline int __cdecl _vsprintf_p_l(
                 char*       const _Buffer,
                                         size_t      const _BufferCount,
        char const* const _Format,
                                     _locale_t   const _Locale,
                                            va_list           _ArgList
    )



{
    int const _Result = __stdio_common_vsprintf_p(
        (*__local_stdio_printf_options()),
        _Buffer, _BufferCount, _Format, _Locale, _ArgList);

    return _Result < 0 ? -1 : _Result;
}


 

__inline int __cdecl _vsprintf_p(
       char*       const _Buffer,
                               size_t      const _BufferCount,
        char const* const _Format,
                                  va_list           _ArgList
    )



{
    return _vsprintf_p_l(_Buffer, _BufferCount, _Format, 0, _ArgList);
}


 

__inline int __cdecl _vsnprintf_s_l(
                 char*       const _Buffer,
                                         size_t      const _BufferCount,
                                         size_t      const _MaxCount,
        char const* const _Format,
                                     _locale_t   const _Locale,
                                            va_list          _ArgList
    )



{
    int const _Result = __stdio_common_vsnprintf_s(
        (*__local_stdio_printf_options()),
        _Buffer, _BufferCount, _MaxCount, _Format, _Locale, _ArgList);

    return _Result < 0 ? -1 : _Result;
}


 

__inline int __cdecl _vsnprintf_s(
       char*       const _Buffer,
                               size_t      const _BufferCount,
                               size_t      const _MaxCount,
        char const* const _Format,
                                  va_list           _ArgList
    )



{
    return _vsnprintf_s_l(_Buffer, _BufferCount, _MaxCount, _Format, 0, _ArgList);
}


extern "C++" { template <size_t _Size> inline   int __cdecl _vsnprintf_s(  char (&_Buffer)[_Size],   size_t _BufferCount,     char const* _Format, va_list _ArgList) throw() { return _vsnprintf_s(_Buffer, _Size, _BufferCount, _Format, _ArgList); } }



     
    
    __inline int __cdecl vsnprintf_s(
           char*       const _Buffer,
                                   size_t      const _BufferCount,
                                   size_t      const _MaxCount,
            char const* const _Format,
                                      va_list           _ArgList
        )



    {
        return _vsnprintf_s_l(_Buffer, _BufferCount, _MaxCount, _Format, 0, _ArgList);
    }

    
    extern "C++" { template <size_t _Size> inline   int __cdecl vsnprintf_s(  char (&_Buffer)[_Size],   size_t _BufferCount,     char const* _Format, va_list _ArgList) throw() { return vsnprintf_s(_Buffer, _Size, _BufferCount, _Format, _ArgList); } }




__inline int __cdecl _vscprintf_l(
        char const* const _Format,
                                     _locale_t   const _Locale,
                                            va_list           _ArgList
    )



{
    int const _Result = __stdio_common_vsprintf(
        (*__local_stdio_printf_options()) | (1ULL << 1),
        0, 0, _Format, _Locale, _ArgList);

    return _Result < 0 ? -1 : _Result;
}


 
__inline int __cdecl _vscprintf(
        char const* const _Format,
                                  va_list           _ArgList
    )



{
    return _vscprintf_l(_Format, 0, _ArgList);
}



__inline int __cdecl _vscprintf_p_l(
        char const* const _Format,
                                     _locale_t   const _Locale,
                                            va_list           _ArgList
    )



{
    int const _Result = __stdio_common_vsprintf_p(
        (*__local_stdio_printf_options()) | (1ULL << 1),
        0, 0, _Format, _Locale, _ArgList);

    return _Result < 0 ? -1 : _Result;
}


 
__inline int __cdecl _vscprintf_p(
        char const* const _Format,
                                  va_list           _ArgList
    )



{
    return _vscprintf_p_l(_Format, 0, _ArgList);
}



__inline int __cdecl _vsnprintf_c_l(
                   char*       const _Buffer,
                                         size_t      const _BufferCount,
        char const* const _Format,
                                     _locale_t   const _Locale,
                                            va_list           _ArgList
    )



{
    int const _Result = __stdio_common_vsprintf(
        (*__local_stdio_printf_options()),
        _Buffer, _BufferCount, _Format, _Locale, _ArgList);

    return _Result < 0 ? -1 : _Result;
}


 

__inline int __cdecl _vsnprintf_c(
         char*       const _Buffer,
                               size_t      const _BufferCount,
        char const* const _Format,
                                  va_list           _ArgList
    )



{
    return _vsnprintf_c_l(_Buffer, _BufferCount, _Format, 0, _ArgList);
}


 __declspec(deprecated("This function or variable may be unsafe. Consider using " "_sprintf_s_l" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__inline int __cdecl _sprintf_l(
                         char*       const _Buffer,
        char const* const _Format,
                                     _locale_t   const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));

    #pragma warning(push)
    #pragma warning(disable: 4996) 
    _Result = _vsprintf_l(_Buffer, _Format, _Locale, _ArgList);
    #pragma warning(pop)

    ((void)(_ArgList = (va_list)0));
    return _Result;
}


 

__inline int __cdecl sprintf(
               char*       const _Buffer,
        char const* const _Format,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));

    #pragma warning(push)
    #pragma warning(disable: 4996) 
    _Result = _vsprintf_l(_Buffer, _Format, 0, _ArgList);
    #pragma warning(pop)

    ((void)(_ArgList = (va_list)0));
    return _Result;
}


#pragma warning(push)
#pragma warning(disable: 4996)
__declspec(deprecated("This function or variable may be unsafe. Consider using " "sprintf_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))   int __cdecl sprintf(  char *_Buffer,  char const* _Format, ...); __declspec(deprecated("This function or variable may be unsafe. Consider using " "vsprintf_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))   int __cdecl vsprintf(  char *_Buffer,  char const* _Format, va_list _Args);
#pragma warning(pop)

 

__inline int __cdecl _sprintf_s_l(
                 char*       const _Buffer,
                                         size_t      const _BufferCount,
        char const* const _Format,
                                     _locale_t   const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vsprintf_s_l(_Buffer, _BufferCount, _Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}




     
    
    __inline int __cdecl sprintf_s(
           char*       const _Buffer,
                                   size_t      const _BufferCount,
            char const* const _Format,
        ...)



    {
        int _Result;
        va_list _ArgList;
        ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
        _Result = _vsprintf_s_l(_Buffer, _BufferCount, _Format, 0, _ArgList);
        ((void)(_ArgList = (va_list)0));
        return _Result;
    }




extern "C++" { __pragma(warning(push)); __pragma(warning(disable: 4793)); template <size_t _Size> inline int __cdecl sprintf_s(  char (&_Buffer)[_Size],     char const* _Format, ...) throw() { va_list _ArgList; ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format))))); return vsprintf_s(_Buffer, _Size, _Format, _ArgList); } __pragma(warning(pop)); }

 

__inline int __cdecl _sprintf_p_l(
                 char*       const _Buffer,
                                         size_t      const _BufferCount,
        char const* const _Format,
                                     _locale_t   const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vsprintf_p_l(_Buffer, _BufferCount, _Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}


 

__inline int __cdecl _sprintf_p(
       char*       const _Buffer,
                               size_t      const _BufferCount,
        char const* const _Format,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
    _Result = _vsprintf_p_l(_Buffer, _BufferCount, _Format, 0, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}


 
 __declspec(deprecated("This function or variable may be unsafe. Consider using " "_snprintf_s_l" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__inline int __cdecl _snprintf_l(
                   char*       const _Buffer,
                                         size_t      const _BufferCount,
        char const* const _Format,
                                     _locale_t   const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));

    #pragma warning(push)
    #pragma warning(disable: 4996) 
    _Result = _vsnprintf_l(_Buffer, _BufferCount, _Format, _Locale, _ArgList);
    #pragma warning(pop)

    ((void)(_ArgList = (va_list)0));
    return _Result;
}













 
 
__inline int __cdecl snprintf(
       char*       const _Buffer,
                               size_t      const _BufferCount,
        char const* const _Format,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
#pragma warning(suppress:28719)    
    _Result = vsnprintf(_Buffer, _BufferCount, _Format, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}


 
__inline int __cdecl _snprintf(
        char*       const _Buffer,
                                          size_t      const _BufferCount,
                   char const* const _Format,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
#pragma warning(suppress:28719)    
    _Result = _vsnprintf(_Buffer, _BufferCount, _Format, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}


__declspec(deprecated("This function or variable may be unsafe. Consider using " "_snprintf_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))    int __cdecl _snprintf(    char *_Buffer,   size_t _BufferCount,     char const* _Format, ...); __declspec(deprecated("This function or variable may be unsafe. Consider using " "_vsnprintf_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))    int __cdecl _vsnprintf(    char *_Buffer,   size_t _BufferCount,     char const* _Format, va_list _Args);

 

__inline int __cdecl _snprintf_c_l(
                   char*       const _Buffer,
                                         size_t      const _BufferCount,
        char const* const _Format,
                                     _locale_t   const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vsnprintf_c_l(_Buffer, _BufferCount, _Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}


 

__inline int __cdecl _snprintf_c(
         char*       const _Buffer,
                               size_t      const _BufferCount,
        char const* const _Format,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
    _Result = _vsnprintf_c_l(_Buffer, _BufferCount, _Format, 0, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}


 

__inline int __cdecl _snprintf_s_l(
                 char*       const _Buffer,
                                         size_t      const _BufferCount,
                                         size_t      const _MaxCount,
        char const* const _Format,
                                     _locale_t   const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vsnprintf_s_l(_Buffer, _BufferCount, _MaxCount, _Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}


 

__inline int __cdecl _snprintf_s(
       char*       const _Buffer,
                               size_t      const _BufferCount,
                               size_t      const _MaxCount,
        char const* const _Format,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
    _Result = _vsnprintf_s_l(_Buffer, _BufferCount, _MaxCount, _Format, 0, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}


extern "C++" { __pragma(warning(push)); __pragma(warning(disable: 4793)); template <size_t _Size> inline   int __cdecl _snprintf_s(  char (&_Buffer)[_Size],   size_t _BufferCount,     char const* _Format, ...) throw() { va_list _ArgList; ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format))))); return _vsnprintf_s(_Buffer, _Size, _BufferCount, _Format, _ArgList); } __pragma(warning(pop)); }


__inline int __cdecl _scprintf_l(
        char const* const _Format,
                                     _locale_t   const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vscprintf_l(_Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}


 
__inline int __cdecl _scprintf(
        char const* const _Format,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
    _Result = _vscprintf_l(_Format, 0, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}



__inline int __cdecl _scprintf_p_l(
        char const* const _Format,
                                     _locale_t   const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vscprintf_p_l(_Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}


 
__inline int __cdecl _scprintf_p(
        char const* const _Format,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
    _Result = _vscprintf_p(_Format, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}







__declspec(dllimport) int __cdecl __stdio_common_vsscanf(
                                        unsigned __int64 _Options,
              char const*      _Buffer,
                                        size_t           _BufferCount,
        char const*      _Format,
                                    _locale_t        _Locale,
                                           va_list          _ArgList
    );


__inline int __cdecl _vsscanf_l(
                             char const* const _Buffer,
        char const* const _Format,
                           _locale_t   const _Locale,
                                  va_list           _ArgList
    )



{
    return __stdio_common_vsscanf(
        (*__local_stdio_scanf_options ()),
        _Buffer, (size_t)-1, _Format, _Locale, _ArgList);
}



__inline int __cdecl vsscanf(
                             char const* const _Buffer,
        char const* const _Format,
                                  va_list           _ArgList
    )



{
    return _vsscanf_l(_Buffer, _Format, 0, _ArgList);
}



__inline int __cdecl _vsscanf_s_l(
                             char const* const _Buffer,
        char const* const _Format,
                           _locale_t   const _Locale,
                                  va_list           _ArgList
    )



{
    return __stdio_common_vsscanf(
        (*__local_stdio_scanf_options ()) | (1ULL << 0),
        _Buffer, (size_t)-1, _Format, _Locale, _ArgList);
}




    #pragma warning(push)
    #pragma warning(disable:6530)

    
    __inline int __cdecl vsscanf_s(
                                 char const* const _Buffer,
            char const* const _Format,
                                      va_list           _ArgList
        )



    {
        return _vsscanf_s_l(_Buffer, _Format, 0, _ArgList);
    }


    extern "C++" { template <size_t _Size> inline int __cdecl vsscanf_s(  char const (&_Buffer)[_Size],     char const* _Format, va_list _ArgList) throw() { return vsscanf_s(_Buffer, _Size, _Format, _ArgList); } }
   
    #pragma warning(pop)



 __declspec(deprecated("This function or variable may be unsafe. Consider using " "_sscanf_s_l" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__inline int __cdecl _sscanf_l(
                                      char const* const _Buffer,
        char const* const _Format,
                                    _locale_t   const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vsscanf_l(_Buffer, _Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}


  __declspec(deprecated("This function or variable may be unsafe. Consider using " "sscanf_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__inline int __cdecl sscanf(
                            char const* const _Buffer,
        char const* const _Format,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));
    _Result = _vsscanf_l(_Buffer, _Format, 0, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}



__inline int __cdecl _sscanf_s_l(
                                        char const* const _Buffer,
        char const* const _Format,
                                      _locale_t   const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));
    _Result = _vsscanf_s_l(_Buffer, _Format, _Locale, _ArgList);
    ((void)(_ArgList = (va_list)0));
    return _Result;
}




    
    __inline int __cdecl sscanf_s(
                                  char const* const _Buffer,
            char const* const _Format,
        ...)



    {
        int _Result;
        va_list _ArgList;
        ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));

        #pragma warning(push)
        #pragma warning(disable: 4996) 
        _Result = vsscanf_s(_Buffer, _Format, _ArgList);
        #pragma warning(pop)

        ((void)(_ArgList = (va_list)0));
        return _Result;
    }




#pragma warning(push)
#pragma warning(disable:6530)

 __declspec(deprecated("This function or variable may be unsafe. Consider using " "_snscanf_s_l" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__inline int __cdecl _snscanf_l(
        char const* const _Buffer,
                                        size_t      const _BufferCount,
        char const* const _Format,
                                    _locale_t   const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));

    _Result = __stdio_common_vsscanf(
        (*__local_stdio_scanf_options ()),
        _Buffer, _BufferCount, _Format, _Locale, _ArgList);

    ((void)(_ArgList = (va_list)0));
    return _Result;
}


 __declspec(deprecated("This function or variable may be unsafe. Consider using " "_snscanf_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__inline int __cdecl _snscanf(
        char const* const _Buffer,
                                        size_t      const _BufferCount,
                  char const* const _Format,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));

    _Result = __stdio_common_vsscanf(
        (*__local_stdio_scanf_options ()),
        _Buffer, _BufferCount, _Format, 0, _ArgList);

    ((void)(_ArgList = (va_list)0));
    return _Result;
}




__inline int __cdecl _snscanf_s_l(
          char const* const _Buffer,
                                          size_t      const _BufferCount,
        char const* const _Format,
                                      _locale_t   const _Locale,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Locale)>(), ((void)(__va_start(&_ArgList, _Locale)))));

    _Result = __stdio_common_vsscanf(
        (*__local_stdio_scanf_options ()) | (1ULL << 0),
        _Buffer, _BufferCount, _Format, _Locale, _ArgList);

    ((void)(_ArgList = (va_list)0));
    return _Result;
}



__inline int __cdecl _snscanf_s(
        char const* const _Buffer,
                                        size_t      const _BufferCount,
                char const* const _Format,
    ...)



{
    int _Result;
    va_list _ArgList;
    ((void)(__vcrt_va_start_verify_argument_type<decltype(_Format)>(), ((void)(__va_start(&_ArgList, _Format)))));

    _Result = __stdio_common_vsscanf(
        (*__local_stdio_scanf_options ()) | (1ULL << 0),
        _Buffer, _BufferCount, _Format, 0, _ArgList);

    ((void)(_ArgList = (va_list)0));
    return _Result;
}


#pragma warning(pop)














    

    




    __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_tempnam" ". See online help for details."))
    __declspec(dllimport) char* __cdecl tempnam(
          char const* _Directory,
          char const* _FilePrefix
        );

    



     __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_fcloseall" ". See online help for details.")) __declspec(dllimport) int   __cdecl fcloseall(void);
          __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_fdopen" ". See online help for details."))    __declspec(dllimport) FILE* __cdecl fdopen(  int _FileHandle,   char const* _Format);
     __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_fgetchar" ". See online help for details."))  __declspec(dllimport) int   __cdecl fgetchar(void);
          __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_fileno" ". See online help for details."))    __declspec(dllimport) int   __cdecl fileno(  FILE* _Stream);
     __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_flushall" ". See online help for details."))  __declspec(dllimport) int   __cdecl flushall(void);
     __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_fputchar" ". See online help for details."))  __declspec(dllimport) int   __cdecl fputchar(  int _Ch);
          __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_getw" ". See online help for details."))      __declspec(dllimport) int   __cdecl getw(  FILE* _Stream);
     __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_putw" ". See online help for details."))      __declspec(dllimport) int   __cdecl putw(  int _Ch,   FILE* _Stream);
          __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_rmtmp" ". See online help for details."))     __declspec(dllimport) int   __cdecl rmtmp(void);





} __pragma(pack(pop))





 
 
 
 
 
 
 

 
 
 
 

  

typedef FILE FILE;

 
namespace std {
using :: FILE; using :: _Mbstatet;

using :: size_t; using :: fpos_t; using :: FILE;
using :: clearerr; using :: fclose; using :: feof;
using :: ferror; using :: fflush; using :: fgetc;
using :: fgetpos; using :: fgets; using :: fopen;
using :: fprintf; using :: fputc; using :: fputs;
using :: fread; using :: freopen; using :: fscanf;
using :: fseek; using :: fsetpos; using :: ftell;
using :: fwrite; using :: getc; using :: getchar;
using :: perror;
using :: putc; using :: putchar;
using :: printf; using :: puts; using :: remove;
using :: rename; using :: rewind; using :: scanf;
using :: setbuf; using :: setvbuf; using :: sprintf;
using :: sscanf; using :: tmpfile; using :: tmpnam;
using :: ungetc; using :: vfprintf; using :: vprintf;
using :: vsprintf;

using :: snprintf; using :: vsnprintf;
using :: vfscanf; using :: vscanf; using :: vsscanf;
}
 










#pragma once










 







#pragma once














#pragma once







__pragma(pack(push, 8)) extern "C" {



 
__declspec(dllimport) int __cdecl _memicmp(
      void const* _Buf1,
      void const* _Buf2,
                             size_t      _Size
    );

 
__declspec(dllimport) int __cdecl _memicmp_l(
      void const* _Buf1,
      void const* _Buf2,
                             size_t      _Size,
                         _locale_t   _Locale
    );





    














    



















    __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_memccpy" ". See online help for details."))
    __declspec(dllimport) void* __cdecl memccpy(
          void*       _Dst,
            void const* _Src,
                                   int         _Val,
                                   size_t      _Size
        );

      __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_memicmp" ". See online help for details."))
    __declspec(dllimport) int __cdecl memicmp(
          void const* _Buf1,
          void const* _Buf2,
                                 size_t      _Size
        );





    extern "C++"  
    inline void* __cdecl memchr(
          void*  _Pv,
                              int    _C,
                              size_t _N
        )
    {
        void const* const _Pvc = _Pv;
        return const_cast<void*>(memchr(_Pvc, _C, _N));
    }




} __pragma(pack(pop))








__pragma(pack(push, 8)) extern "C" {







     
    __declspec(dllimport) errno_t __cdecl strcpy_s(
          char*       _Destination,
                                  rsize_t     _SizeInBytes,
                                char const* _Source
        );

    
    __declspec(dllimport) errno_t __cdecl strcat_s(
          char*       _Destination,
                                     rsize_t     _SizeInBytes,
                                   char const* _Source
        );

    
    __declspec(dllimport) errno_t __cdecl strerror_s(
          char*  _Buffer,
                                  size_t _SizeInBytes,
                                  int    _ErrorNumber);

    
    __declspec(dllimport) errno_t __cdecl strncat_s(
          char*       _Destination,
                                     rsize_t     _SizeInBytes,
               char const* _Source,
                                     rsize_t     _MaxCount
        );

    
    __declspec(dllimport) errno_t __cdecl strncpy_s(
          char*       _Destination,
                                  rsize_t     _SizeInBytes,
            char const* _Source,
                                  rsize_t     _MaxCount
        );

     
    __declspec(dllimport) char*  __cdecl strtok_s(
                          char*       _String,
                                 char const* _Delimiter,
            char**      _Context
        );



__declspec(dllimport) void* __cdecl _memccpy(
      void*       _Dst,
                                   void const* _Src,
                                   int         _Val,
                                   size_t      _MaxCount
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl strcat_s(char (&_Destination)[_Size],   char const* _Source) throw() { return strcat_s(_Destination, _Size, _Source); } }



    __declspec(deprecated("This function or variable may be unsafe. Consider using " "strcat_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))  char* __cdecl strcat( char *_Destination,  char const* _Source);



 
int __cdecl strcmp(
      char const* _Str1,
      char const* _Str2
    );

 
__declspec(dllimport) int __cdecl _strcmpi(
      char const* _String1,
      char const* _String2
    );

 
__declspec(dllimport) int __cdecl strcoll(
      char const* _String1,
      char const* _String2
    );

 
__declspec(dllimport) int __cdecl _strcoll_l(
        char const* _String1,
        char const* _String2,
      _locale_t   _Locale
    );

char* __cdecl strcpy(
      char*       _Dest,
                                            char const* _Source
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl strcpy_s(  char (&_Destination)[_Size],   char const* _Source) throw() { return strcpy_s(_Destination, _Size, _Source); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "strcpy_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))  char* __cdecl strcpy( char *_Destination,  char const* _Source);

 
__declspec(dllimport) size_t __cdecl strcspn(
      char const* _Str,
      char const* _Control
    );






 
__declspec(dllimport) __declspec(allocator) char* __cdecl _strdup(
      char const* _Source
    );





 
 
  __declspec(deprecated("This function or variable may be unsafe. Consider using " "_strerror_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__declspec(dllimport) char*  __cdecl _strerror(
      char const* _ErrorMessage
    );


__declspec(dllimport) errno_t __cdecl _strerror_s(
      char*       _Buffer,
                              size_t      _SizeInBytes,
                        char const* _ErrorMessage
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl _strerror_s(char (&_Buffer)[_Size],   char const* _ErrorMessage) throw() { return _strerror_s(_Buffer, _Size, _ErrorMessage); } }

 
  __declspec(deprecated("This function or variable may be unsafe. Consider using " "strerror_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__declspec(dllimport) char* __cdecl strerror(
      int _ErrorMessage
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl strerror_s(char (&_Buffer)[_Size],   int _ErrorMessage) throw() { return strerror_s(_Buffer, _Size, _ErrorMessage); } }

 
__declspec(dllimport) int __cdecl _stricmp(
      char const* _String1,
      char const* _String2
    );

 
__declspec(dllimport) int __cdecl _stricoll(
      char const* _String1,
      char const* _String2
    );

 
__declspec(dllimport) int __cdecl _stricoll_l(
        char const* _String1,
        char const* _String2,
      _locale_t   _Locale
    );

 
__declspec(dllimport) int __cdecl _stricmp_l(
        char const* _String1,
        char const* _String2,
      _locale_t   _Locale
    );

 
size_t __cdecl strlen(
      char const* _Str
    );


__declspec(dllimport) errno_t __cdecl _strlwr_s(
      char*  _String,
                          size_t _Size
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl _strlwr_s(  char (&_String)[_Size]) throw() { return _strlwr_s(_String, _Size); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "_strlwr_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport) char* __cdecl _strlwr( char *_String);


__declspec(dllimport) errno_t __cdecl _strlwr_s_l(
      char*     _String,
                          size_t    _Size,
                      _locale_t _Locale
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl _strlwr_s_l(  char (&_String)[_Size],   _locale_t _Locale) throw() { return _strlwr_s_l(_String, _Size, _Locale); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "_strlwr_s_l" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport) char* __cdecl _strlwr_l(  char *_String,   _locale_t _Locale);

__declspec(dllimport) char* __cdecl strncat(
      char*       _Dest,
        char const* _Source,
                           size_t      _Count
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl strncat_s(  char (&_Destination)[_Size],   char const* _Source,   size_t _Count) throw() { return strncat_s(_Destination, _Size, _Source, _Count); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "strncat_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport) char* __cdecl strncat(  char *_Destination,   char const* _Source,   size_t _Count);

 
__declspec(dllimport) int __cdecl strncmp(
      char const* _Str1,
      char const* _Str2,
                            size_t      _MaxCount
    );

 
__declspec(dllimport) int __cdecl _strnicmp(
      char const* _String1,
      char const* _String2,
                            size_t      _MaxCount
    );

 
__declspec(dllimport) int __cdecl _strnicmp_l(
      char const* _String1,
      char const* _String2,
                            size_t      _MaxCount,
                        _locale_t   _Locale
    );

 
__declspec(dllimport) int __cdecl _strnicoll(
      char const* _String1,
      char const* _String2,
                            size_t      _MaxCount
    );

 
__declspec(dllimport) int __cdecl _strnicoll_l(
      char const* _String1,
      char const* _String2,
                            size_t      _MaxCount,
                        _locale_t   _Locale
    );

 
__declspec(dllimport) int __cdecl _strncoll(
      char const* _String1,
      char const* _String2,
                            size_t      _MaxCount
    );

 
__declspec(dllimport) int __cdecl _strncoll_l(
      char const* _String1,
      char const* _String2,
                            size_t      _MaxCount,
                        _locale_t   _Locale
    );

__declspec(dllimport) size_t __cdecl __strncnt(
      char const* _String,
                         size_t      _Count
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl strncpy_s(char (&_Destination)[_Size],   char const* _Source,   size_t _Count) throw() { return strncpy_s(_Destination, _Size, _Source, _Count); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "strncpy_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport) char* __cdecl strncpy(    char *_Destination,   char const* _Source,   size_t _Count);

 


__declspec(dllimport) size_t __cdecl strnlen(
      char const* _String,
                            size_t      _MaxCount
    );



     
    
    
    static __inline size_t __cdecl strnlen_s(
          char const* _String,
                                size_t      _MaxCount
        )
    {
        return _String == 0 ? 0 : strnlen(_String, _MaxCount);
    }



__declspec(dllimport) char* __cdecl _strnset(
      char*  _Dest,
                           int    _Val,
                           size_t _Count
    );


__declspec(dllimport) errno_t __cdecl _strnset_s(
      char*  _String,
                                 size_t _SizeInBytes,
                                 int    _Value,
                                 size_t _MaxCount
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl _strnset_s(  char (&_Destination)[_Size],   int _Value,   size_t _Count) throw() { return _strnset_s(_Destination, _Size, _Value, _Count); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "_strnset_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport) char* __cdecl _strnset(  char *_Destination,   int _Value,   size_t _Count);

 
__declspec(dllimport) char const* __cdecl strpbrk(
      char const* _Str,
      char const* _Control
    );

__declspec(dllimport) char* __cdecl _strrev(
      char* _Str
    );


__declspec(dllimport) errno_t __cdecl _strset_s(
      char*  _Destination,
                                     size_t _DestinationSize,
                                     int    _Value
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl _strset_s(  char (&_Destination)[_Size],   int _Value) throw() { return _strset_s(_Destination, _Size, _Value); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "_strset_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))  char* __cdecl _strset( char *_Destination,  int _Value);

char* __cdecl _strset(
      char* _Dest,
           int   _Value
    );

 
__declspec(dllimport) size_t __cdecl strspn(
      char const* _Str,
      char const* _Control
    );

  __declspec(deprecated("This function or variable may be unsafe. Consider using " "strtok_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details."))
__declspec(dllimport) char* __cdecl strtok(
      char*       _String,
             char const* _Delimiter
    );


__declspec(dllimport) errno_t __cdecl _strupr_s(
      char*  _String,
                          size_t _Size
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl _strupr_s(  char (&_String)[_Size]) throw() { return _strupr_s(_String, _Size); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "_strupr_s" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport) char* __cdecl _strupr( char *_String);


__declspec(dllimport) errno_t __cdecl _strupr_s_l(
      char*     _String,
                          size_t    _Size,
                      _locale_t _Locale
    );

extern "C++" { template <size_t _Size> inline errno_t __cdecl _strupr_s_l(  char (&_String)[_Size],   _locale_t _Locale) throw() { return _strupr_s_l(_String, _Size, _Locale); } }

__declspec(deprecated("This function or variable may be unsafe. Consider using " "_strupr_s_l" " instead. To disable deprecation, use _CRT_SECURE_NO_WARNINGS. " "See online help for details.")) __declspec(dllimport) char* __cdecl _strupr_l(  char *_String,   _locale_t _Locale);

 

__declspec(dllimport) size_t __cdecl strxfrm(
        char*       _Destination,
                                         char const* _Source,
                 size_t      _MaxCount
    );

 

__declspec(dllimport) size_t __cdecl _strxfrm_l(
        char*       _Destination,
                                         char const* _Source,
                 size_t      _MaxCount,
                                       _locale_t   _Locale
    );




extern "C++"
{
     
    inline char* __cdecl strchr(  char* const _String,   int const _Ch)
    {
        return const_cast<char*>(strchr(static_cast<char const*>(_String), _Ch));
    }

     
    inline char* __cdecl strpbrk(  char* const _String,   char const* const _Control)
    {
        return const_cast<char*>(strpbrk(static_cast<char const*>(_String), _Control));
    }

     
    inline char* __cdecl strrchr(  char* const _String,   int const _Ch)
    {
        return const_cast<char*>(strrchr(static_cast<char const*>(_String), _Ch));
    }

       
    inline char* __cdecl strstr(  char* const _String,   char const* const _SubString)
    {
        return const_cast<char*>(strstr(static_cast<char const*>(_String), _SubString));
    }
}






    




      __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_strdup" ". See online help for details."))
    __declspec(dllimport) char* __cdecl strdup(
          char const* _String
        );

    



    
      __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_strcmpi" ". See online help for details."))
    __declspec(dllimport) int __cdecl strcmpi(
          char const* _String1,
          char const* _String2
        );

      __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_stricmp" ". See online help for details."))
    __declspec(dllimport) int __cdecl stricmp(
          char const* _String1,
          char const* _String2
        );

    __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_strlwr" ". See online help for details."))
    __declspec(dllimport) char* __cdecl strlwr(
          char* _String
        );

      __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_strnicmp" ". See online help for details."))
    __declspec(dllimport) int __cdecl strnicmp(
          char const* _String1,
          char const* _String2,
                                size_t      _MaxCount
        );

    __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_strnset" ". See online help for details."))
    __declspec(dllimport) char* __cdecl strnset(
          char*  _String,
                                  int    _Value,
                                  size_t _MaxCount
        );

    __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_strrev" ". See online help for details."))
    __declspec(dllimport) char* __cdecl strrev(
          char* _String
        );

    __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_strset" ". See online help for details."))
    char* __cdecl strset(
          char* _String,
               int   _Value);

    __declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C " "and C++ conformant name: " "_strupr" ". See online help for details."))
    __declspec(dllimport) char* __cdecl strupr(
          char* _String
        );





} __pragma(pack(pop))





 
namespace std {
using :: size_t; using :: memchr; using :: memcmp;
using :: memcpy; using :: memmove; using :: memset;
using :: strcat; using :: strchr; using :: strcmp;
using :: strcoll; using :: strcpy; using :: strcspn;
using :: strerror; using :: strlen; using :: strncat;
using :: strncmp; using :: strncpy; using :: strpbrk;
using :: strrchr; using :: strspn; using :: strstr;
using :: strtok; using :: strxfrm;
}
 



















#pragma once











#pragma once




extern "C++" {

#pragma pack(push, 8)






         
    __declspec(allocator) void* __cdecl operator new(
            size_t      _Size,
            int         _BlockUse,
          char const* _FileName,
            int         _LineNumber
        );

         
    __declspec(allocator) void* __cdecl operator new[](
            size_t      _Size,
            int         _BlockUse,
          char const* _FileName,
            int         _LineNumber
        );

    void __cdecl operator delete(
        void*       _Block,
        int         _BlockUse,
        char const* _FileName,
        int         _LineNumber
        ) throw();

    void __cdecl operator delete[](
        void*       _Block,
        int         _BlockUse,
        char const* _FileName,
        int         _LineNumber
        ) throw();





#pragma pack(pop)

} 



__pragma(pack(push, 8)) extern "C" {



typedef void* _HFILE; 

























typedef int (__cdecl* _CRT_REPORT_HOOK )(int, char*,    int*);
typedef int (__cdecl* _CRT_REPORT_HOOKW)(int, wchar_t*, int*);





typedef int (__cdecl* _CRT_ALLOC_HOOK)(int, void*, size_t, int, long, unsigned char const*, int);























































typedef void (__cdecl* _CRT_DUMP_CLIENT)(void*, size_t);





struct _CrtMemBlockHeader;

typedef struct _CrtMemState
{
    struct _CrtMemBlockHeader* pBlockHeader;
    size_t lCounts[5];
    size_t lSizes[5];
    size_t lHighWaterCount;
    size_t lTotalCount;
} _CrtMemState;



    
    

    
    

    
    
    
    
    
    
    
    
    
    
    
    
    



































































































    
    
    
    
    
    
    

    
    
    
    
    
    
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    





























































































































































































































































































































    
    
    
    
    
    
    





























































































    

    
        
    

    
        
    

    
        
    

    
    

    
    

    
    

    
    


























































    




























} __pragma(pack(pop))



 #pragma pack(push,8)
 #pragma warning(push,3)
 
 

namespace std {
		

typedef _Longlong streamoff;
typedef _Longlong streamsize;

  
  

  



extern __declspec(dllimport)  const streamoff _BADOFF;
  

		
template<class _Statetype>
	class fpos
	{	
	typedef fpos<_Statetype> _Myt;

public:
	 fpos(streamoff _Off = 0)
		: _Myoff(_Off), _Fpos(0), _Mystate()
		{	
		}

	 fpos(_Statetype _State, fpos_t _Fileposition)
		: _Myoff(0), _Fpos(_Fileposition), _Mystate(_State)
		{	
		}

	_Statetype  state() const
		{	
		return (_Mystate);
		}

	void  state(_Statetype _State)
		{	
		_Mystate = _State;
		}

	fpos_t  seekpos() const
		{	
		return (_Fpos);
		}

	 operator streamoff() const
		{	
		return ((streamoff)(_Myoff + ((long long)(_Fpos))));
		}

	streamoff  operator-(const _Myt& _Right) const
		{	
		return ((streamoff)*this - (streamoff)_Right);
		}

	_Myt&  operator+=(streamoff _Off)
		{	
		_Myoff += _Off;
		return (*this);
		}

	_Myt&  operator-=(streamoff _Off)
		{	
		_Myoff -= _Off;
		return (*this);
		}

	_Myt  operator+(streamoff _Off) const
		{	
		_Myt _Tmp = *this;
		return (_Tmp += _Off);
		}

	_Myt  operator-(streamoff _Off) const
		{	
		_Myt _Tmp = *this;
		return (_Tmp -= _Off);
		}

	bool  operator==(const _Myt& _Right) const
		{	
		return ((streamoff)*this == (streamoff)_Right);
		}

	bool  operator==(streamoff _Right) const
		{	
		return ((streamoff)*this == _Right);
		}

	bool  operator!=(const _Myt& _Right) const
		{	
		return (!(*this == _Right));
		}

private:
	streamoff _Myoff;	
	fpos_t _Fpos;	
	_Statetype _Mystate;	
	};

 

 
 

typedef fpos<_Mbstatet> streampos;

typedef streampos wstreampos;

		
template<class _Elem,
	class _Int_type>
	struct _Char_traits
	{	
	typedef _Elem char_type;
	typedef _Int_type int_type;
	typedef streampos pos_type;
	typedef streamoff off_type;
	typedef _Mbstatet state_type;

	static int __cdecl compare(
		  const _Elem *_First1,
		  const _Elem *_First2, size_t _Count)
		{	
		for (; 0 < _Count; --_Count, ++_First1, ++_First2)
			if (!eq(*_First1, *_First2))
				return (lt(*_First1, *_First2) ? -1 : +1);
		return (0);
		}

	static size_t __cdecl length(  const _Elem *_First)
		{	
		size_t _Count;
		for (_Count = 0; !eq(*_First, _Elem()); ++_First)
			++_Count;
		return (_Count);
		}

	static _Elem *__cdecl copy(
		  _Elem *_First1,
		  const _Elem *_First2, size_t _Count)
		{	
		_Elem *_Next = _First1;
		for (; 0 < _Count; --_Count, ++_Next, ++_First2)
			assign(*_Next, *_First2);
		return (_First1);
		}

	static _Elem *__cdecl _Copy_s(
		  _Elem *_First1, size_t _Dest_size,
		  const _Elem *_First2, size_t _Count)
		{	
		{ if (!(_Count <= _Dest_size)) { ((void)0); ::_invalid_parameter_noinfo_noreturn(); return (0); } };
		return (copy(_First1, _First2, _Count));
		}

	static const _Elem *__cdecl find(
		  const _Elem *_First,
		size_t _Count, const _Elem& _Ch)
		{	
		for (; 0 < _Count; --_Count, ++_First)
			if (eq(*_First, _Ch))
				return (_First);
		return (0);
		}

	static _Elem *__cdecl move(
		  _Elem *_First1,
		  const _Elem *_First2, size_t _Count)
		{	
		_Elem *_Next = _First1;
		if (_First2 < _Next && _Next < _First2 + _Count)
			for (_Next += _Count, _First2 += _Count; 0 < _Count; --_Count)
				assign(*--_Next, *--_First2);
		else
			for (; 0 < _Count; --_Count, ++_Next, ++_First2)
				assign(*_Next, *_First2);
		return (_First1);
		}

	static _Elem *__cdecl assign(
		  _Elem *_First,
		size_t _Count, _Elem _Ch)
		{	
		_Elem *_Next = _First;
		for (; 0 < _Count; --_Count, ++_Next)
			assign(*_Next, _Ch);
		return (_First);
		}

	static void __cdecl assign(_Elem& _Left, const _Elem& _Right) noexcept
		{	
		_Left = _Right;
		}

	static constexpr bool __cdecl eq(const _Elem& _Left,
		const _Elem& _Right) noexcept
		{	
		return (_Left == _Right);
		}

	static constexpr bool __cdecl lt(const _Elem& _Left,
		const _Elem& _Right) noexcept
		{	
		return (_Left < _Right);
		}

	static constexpr _Elem __cdecl to_char_type(
		const int_type& _Meta) noexcept
		{	
		return ((_Elem)_Meta);
		}

	static constexpr int_type __cdecl to_int_type(
		const _Elem& _Ch) noexcept
		{	
		return ((int_type)_Ch);
		}

	static constexpr bool __cdecl eq_int_type(const int_type& _Left,
		const int_type& _Right) noexcept
		{	
		return (_Left == _Right);
		}

	static constexpr int_type __cdecl not_eof(
		const int_type& _Meta) noexcept
		{	
		return (_Meta != eof() ? (int_type)_Meta : (int_type)!eof());
		}

	static constexpr int_type __cdecl eof() noexcept
		{	
		return ((int_type)(-1));
		}
	};

		
template<class _Elem>
	struct char_traits
		: public _Char_traits<_Elem, long>
	{	
	};

		
template<>
	struct char_traits<char16_t>
	: public _Char_traits<char16_t, unsigned short>
	{	
	};

typedef streampos u16streampos;

		
template<>
	struct char_traits<char32_t>
	: public _Char_traits<char32_t, unsigned int>
	{	
	};

typedef streampos u32streampos;

		
template<>
	struct char_traits<wchar_t>
	{	
	typedef wchar_t _Elem;
	typedef _Elem char_type;	
	typedef wint_t int_type;
	typedef streampos pos_type;
	typedef streamoff off_type;
	typedef _Mbstatet state_type;

	static int __cdecl compare(const _Elem *_First1, const _Elem *_First2,
		size_t _Count)
		{	
		return (_Count == 0 ? 0
			: :: wmemcmp(_First1, _First2, _Count));
		}

	static size_t __cdecl length(const _Elem *_First)
		{	
		return (*_First == 0 ? 0
			: :: wcslen(_First));
		}

	static _Elem *__cdecl copy(_Elem *_First1, const _Elem *_First2,
		size_t _Count)
		{	
		return (_Count == 0 ? _First1
			: (_Elem *):: wmemcpy(_First1, _First2, _Count));
		}

	static _Elem *__cdecl _Copy_s(
		  _Elem *_First1, size_t _Size_in_words,
		  const _Elem *_First2, size_t _Count)
		{	
		if (0 < _Count)
			::wmemcpy_s((_First1), (_Size_in_words), (_First2), (_Count));
		return (_First1);
		}

	static const _Elem *__cdecl find(const _Elem *_First, size_t _Count,
		const _Elem& _Ch)
		{	
		return (_Count == 0 ? (const _Elem *)0
			: (const _Elem *):: wmemchr(_First, _Ch, _Count));
		}

	static _Elem *__cdecl move(_Elem *_First1, const _Elem *_First2,
		size_t _Count)
		{	
		return (_Count == 0 ? _First1
			: (_Elem *):: wmemmove(_First1, _First2, _Count));
		}

	static _Elem *__cdecl assign(_Elem *_First, size_t _Count,
		_Elem _Ch)
		{	
		return ((_Elem *):: wmemset(_First, _Ch, _Count));
		}

	static void __cdecl assign(_Elem& _Left, const _Elem& _Right) noexcept
		{	
		_Left = _Right;
		}

	static constexpr bool __cdecl eq(const _Elem& _Left,
		const _Elem& _Right) noexcept
		{	
		return (_Left == _Right);
		}

	static constexpr bool __cdecl lt(const _Elem& _Left,
		const _Elem& _Right) noexcept
		{	
		return (_Left < _Right);
		}

	static constexpr _Elem __cdecl to_char_type(
		const int_type& _Meta) noexcept
		{	
		return (_Meta);
		}

	static constexpr int_type __cdecl to_int_type(
		const _Elem& _Ch) noexcept
		{	
		return (_Ch);
		}

	static constexpr bool __cdecl eq_int_type(const int_type& _Left,
		const int_type& _Right) noexcept
		{	
		return (_Left == _Right);
		}

	static constexpr int_type __cdecl not_eof(
		const int_type& _Meta) noexcept
		{	
		return (_Meta != eof() ? _Meta : !eof());
		}

	static constexpr int_type __cdecl eof() noexcept
		{	
		return (((wint_t)(0xFFFF)));
		}
	};

 
		
template<>
	struct char_traits<unsigned short>
	{	
	typedef unsigned short _Elem;
	typedef _Elem char_type;	
	typedef wint_t int_type;
	typedef streampos pos_type;
	typedef streamoff off_type;
	typedef _Mbstatet state_type;

	static int __cdecl compare(const _Elem *_First1, const _Elem *_First2,
		size_t _Count)
		{	
		return (_Count == 0 ? 0
			: :: wmemcmp((const wchar_t *)_First1,
				(const wchar_t *)_First2, _Count));
		}

	static size_t __cdecl length(const _Elem *_First)
		{	
		return (*_First == 0 ? 0
			: :: wcslen((const wchar_t *)_First));
		}

	static _Elem *__cdecl copy(_Elem *_First1, const _Elem *_First2,
		size_t _Count)
		{	
		return (_Count == 0 ? _First1
			: (_Elem *):: wmemcpy((wchar_t *)_First1,
				(const wchar_t *)_First2, _Count));
		}

	static _Elem *__cdecl _Copy_s(
		  _Elem *_First1, size_t _Size_in_words,
		  const _Elem *_First2, size_t _Count)
		{	
		if (0 < _Count)
			::wmemcpy_s(((wchar_t *)_First1), (_Size_in_words), ((const wchar_t *)_First2), (_Count));
		return (_First1);
		}

	static const _Elem *__cdecl find(const _Elem *_First, size_t _Count,
		const _Elem& _Ch)
		{	
		return (_Count == 0 ? (const _Elem *)0
			: (const _Elem *):: wmemchr((const wchar_t *)_First,
				_Ch, _Count));
		}

	static _Elem *__cdecl move(_Elem *_First1, const _Elem *_First2,
		size_t _Count)
		{	
		return (_Count == 0 ? _First1
			: (_Elem *):: wmemmove((wchar_t *)_First1,
				(const wchar_t *)_First2, _Count));
		}

	static _Elem *__cdecl assign(_Elem *_First, size_t _Count,
		_Elem _Ch)
		{	
		return ((_Elem *):: wmemset((wchar_t *)_First, _Ch, _Count));
		}

	static void __cdecl assign(_Elem& _Left, const _Elem& _Right) noexcept
		{	
		_Left = _Right;
		}

	static constexpr bool __cdecl eq(const _Elem& _Left,
		const _Elem& _Right) noexcept
		{	
		return (_Left == _Right);
		}

	static constexpr bool __cdecl lt(const _Elem& _Left,
		const _Elem& _Right) noexcept
		{	
		return (_Left < _Right);
		}

	static constexpr _Elem __cdecl to_char_type(const int_type& _Meta)
		noexcept
		{	
		return (_Meta);
		}

	static constexpr int_type __cdecl to_int_type(const _Elem& _Ch)
		noexcept
		{	
		return (_Ch);
		}

	static constexpr bool __cdecl eq_int_type(const int_type& _Left,
		const int_type& _Right) noexcept
		{	
		return (_Left == _Right);
		}

	static constexpr int_type __cdecl not_eof(const int_type& _Meta)
		noexcept
		{	
		return (_Meta != eof() ? _Meta : !eof());
		}

	static constexpr int_type __cdecl eof() noexcept
		{	
		return (((wint_t)(0xFFFF)));
		}
	};
 

		
template<> struct char_traits<char>
	{	
	typedef char _Elem;
	typedef _Elem char_type;
	typedef int int_type;
	typedef streampos pos_type;
	typedef streamoff off_type;
	typedef _Mbstatet state_type;

	static int __cdecl compare(const _Elem *_First1, const _Elem *_First2,
		size_t _Count)
		{	
		return (_Count == 0 ? 0
			: :: memcmp(_First1, _First2, _Count));
		}

	static size_t __cdecl length(const _Elem *_First)
		{	
		return (*_First == 0 ? 0
			: :: strlen(_First));
		}

	static _Elem *__cdecl copy(_Elem *_First1, const _Elem *_First2,
		size_t _Count)
		{	
		return (_Count == 0 ? _First1
			: (_Elem *):: memcpy(_First1, _First2, _Count));
		}

	static _Elem *__cdecl _Copy_s(
		  _Elem *_First1, size_t _Size_in_bytes,
		  const _Elem *_First2, size_t _Count)
		{	
		if (0 < _Count)
			::memcpy_s((_First1), (_Size_in_bytes), (_First2), (_Count));
		return (_First1);
		}

	static const _Elem *__cdecl find(const _Elem *_First, size_t _Count,
		const _Elem& _Ch)
		{	
		return (_Count == 0 ? (const _Elem *)0
			: (const _Elem *):: memchr(_First, _Ch, _Count));
		}

	static _Elem *__cdecl move(_Elem *_First1, const _Elem *_First2,
		size_t _Count)
		{	
		return (_Count == 0 ? _First1
			: (_Elem *):: memmove(_First1, _First2, _Count));
		}

	static _Elem *__cdecl assign(_Elem *_First, size_t _Count,
		_Elem _Ch)
		{	
		return ((_Elem *):: memset(_First, _Ch, _Count));
		}

	static void __cdecl assign(_Elem& _Left, const _Elem& _Right) noexcept
		{	
		_Left = _Right;
		}

	static constexpr bool __cdecl eq(const _Elem& _Left,
		const _Elem& _Right) noexcept
		{	
		return (_Left == _Right);
		}

	static constexpr bool __cdecl lt(const _Elem& _Left,
		const _Elem& _Right) noexcept
		{	
		return ((unsigned char)_Left < (unsigned char)_Right);
		}

	static constexpr _Elem __cdecl to_char_type(
		const int_type& _Meta) noexcept
		{	
		return ((_Elem)_Meta);
		}

	static constexpr int_type __cdecl to_int_type(
		const _Elem& _Ch) noexcept
		{	
		return ((unsigned char)_Ch);
		}

	static constexpr bool __cdecl eq_int_type(const int_type& _Left,
		const int_type& _Right) noexcept
		{	
		return (_Left == _Right);
		}

	static constexpr int_type __cdecl not_eof(
		const int_type& _Meta) noexcept
		{	
		return (_Meta != eof() ? _Meta : !eof());
		}

	static constexpr int_type __cdecl eof() noexcept
		{	
		return ((-1));
		}
	};

		
template<class _Ty>
	class allocator;
class ios_base;
template<class _Elem,
	class _Traits = char_traits<_Elem> >
	class basic_ios;
template<class _Elem,
	class _Traits = char_traits<_Elem> >
	class istreambuf_iterator;
template<class _Elem,
	class _Traits = char_traits<_Elem> >
	class ostreambuf_iterator;
template<class _Elem,
	class _Traits = char_traits<_Elem> >
	class basic_streambuf;
template<class _Elem,
	class _Traits = char_traits<_Elem> >
	class basic_istream;
template<class _Elem,
	class _Traits = char_traits<_Elem> >
	class basic_ostream;
template<class _Elem,
	class _Traits = char_traits<_Elem> >
	class basic_iostream;
template<class _Elem,
	class _Traits = char_traits<_Elem>,
	class _Alloc = allocator<_Elem> >
	class basic_stringbuf;
template<class _Elem,
	class _Traits = char_traits<_Elem>,
	class _Alloc = allocator<_Elem> >
	class basic_istringstream;
template<class _Elem,
	class _Traits = char_traits<_Elem>,
	class _Alloc = allocator<_Elem> >
	class basic_ostringstream;
template<class _Elem,
	class _Traits = char_traits<_Elem>,
	class _Alloc = allocator<_Elem> >
	class basic_stringstream;
template<class _Elem,
	class _Traits = char_traits<_Elem> >
	class basic_filebuf;
template<class _Elem,
	class _Traits = char_traits<_Elem> >
	class basic_ifstream;
template<class _Elem,
	class _Traits = char_traits<_Elem> >
	class basic_ofstream;
template<class _Elem,
	class _Traits = char_traits<_Elem> >
	class basic_fstream;

 
template<class _Elem,
	class _InIt >
	class num_get;
template<class _Elem,
	class _OutIt >
	class num_put;
template<class _Elem>
	class collate;
 

		
typedef basic_ios<char, char_traits<char> > ios;
typedef basic_streambuf<char, char_traits<char> > streambuf;
typedef basic_istream<char, char_traits<char> > istream;
typedef basic_ostream<char, char_traits<char> > ostream;
typedef basic_iostream<char, char_traits<char> > iostream;
typedef basic_stringbuf<char, char_traits<char>,
	allocator<char> > stringbuf;
typedef basic_istringstream<char, char_traits<char>,
	allocator<char> > istringstream;
typedef basic_ostringstream<char, char_traits<char>,
	allocator<char> > ostringstream;
typedef basic_stringstream<char, char_traits<char>,
	allocator<char> > stringstream;
typedef basic_filebuf<char, char_traits<char> > filebuf;
typedef basic_ifstream<char, char_traits<char> > ifstream;
typedef basic_ofstream<char, char_traits<char> > ofstream;
typedef basic_fstream<char, char_traits<char> > fstream;

		
typedef basic_ios<wchar_t, char_traits<wchar_t> > wios;
typedef basic_streambuf<wchar_t, char_traits<wchar_t> >
	wstreambuf;
typedef basic_istream<wchar_t, char_traits<wchar_t> > wistream;
typedef basic_ostream<wchar_t, char_traits<wchar_t> > wostream;
typedef basic_iostream<wchar_t, char_traits<wchar_t> > wiostream;
typedef basic_stringbuf<wchar_t, char_traits<wchar_t>,
	allocator<wchar_t> > wstringbuf;
typedef basic_istringstream<wchar_t, char_traits<wchar_t>,
	allocator<wchar_t> > wistringstream;
typedef basic_ostringstream<wchar_t, char_traits<wchar_t>,
	allocator<wchar_t> > wostringstream;
typedef basic_stringstream<wchar_t, char_traits<wchar_t>,
	allocator<wchar_t> > wstringstream;
typedef basic_filebuf<wchar_t, char_traits<wchar_t> > wfilebuf;
typedef basic_ifstream<wchar_t, char_traits<wchar_t> > wifstream;
typedef basic_ofstream<wchar_t, char_traits<wchar_t> > wofstream;
typedef basic_fstream<wchar_t, char_traits<wchar_t> > wfstream;

 





















 
typedef num_get<char, istreambuf_iterator<char, char_traits<char> > >
	numget;
typedef num_get<wchar_t, istreambuf_iterator<wchar_t, char_traits<wchar_t> > >
	wnumget;
typedef num_put<char, ostreambuf_iterator<char, char_traits<char> > >
	numput;
typedef num_put<wchar_t, ostreambuf_iterator<wchar_t, char_traits<wchar_t> > >
	wnumput;
typedef collate<char> ncollate;
typedef collate<wchar_t> wcollate;
 
}

 
 #pragma warning(pop)
 #pragma pack(pop)










 #pragma pack(push,8)
 #pragma warning(push,3)
 
 

 #pragma warning(disable: 4180 4512)

namespace std {
		
template<class _FwdIt1,
	class _FwdIt2> inline
	void iter_swap(_FwdIt1 _Left, _FwdIt2 _Right)
	{	
	swap(*_Left, *_Right);
	}

		
template<class _Ty,
	size_t _Size,
	class> inline
	void swap(_Ty (&_Left)[_Size], _Ty (&_Right)[_Size])
		noexcept(_Is_nothrow_swappable<_Ty>::value)
	{	
	if (&_Left != &_Right)
		{	
		_Ty *_First1 = _Left;
		_Ty *_Last1 = _First1 + _Size;
		_Ty *_First2 = _Right;
		for (; _First1 != _Last1; ++_First1, ++_First2)
			::std:: iter_swap(_First1, _First2);
		}
	}

template<class _Ty,
	class> inline
	void swap(_Ty& _Left, _Ty& _Right)
		noexcept(is_nothrow_move_constructible<_Ty>::value && is_nothrow_move_assignable<_Ty>::value)
	{	
	_Ty _Tmp = ::std:: move(_Left);
	_Left = ::std:: move(_Right);
	_Right = ::std:: move(_Tmp);
	}

		
template<class _Ty> inline
	void _Swap_adl(_Ty& _Left, _Ty& _Right)
		noexcept(_Is_nothrow_swappable<_Ty>::value)
	{	
	swap(_Left, _Right);
	}

		
struct piecewise_construct_t
	{	
	};

constexpr piecewise_construct_t piecewise_construct{};

		
template<class...>
	class tuple;

template<class _Ty1,
	class _Ty2>
	struct pair
	{	
	typedef pair<_Ty1, _Ty2> _Myt;
	typedef _Ty1 first_type;
	typedef _Ty2 second_type;

	template<class _Uty1 = _Ty1,
		class _Uty2 = _Ty2,
		class = enable_if_t<is_default_constructible<_Uty1>::value
						&& is_default_constructible<_Uty2>::value> >
		constexpr pair()
		: first(), second()
		{	
		}

	template<class _Uty1 = _Ty1,
		class _Uty2 = _Ty2,
		class = enable_if_t<is_copy_constructible<_Uty1>::value
						&& is_copy_constructible<_Uty2>::value>,
		enable_if_t<is_convertible<const _Uty1&, _Uty1>::value
				&& is_convertible<const _Uty2&, _Uty2>::value, int> = 0>
		constexpr pair(const _Ty1& _Val1, const _Ty2& _Val2)
		: first(_Val1), second(_Val2)
		{	
		}

	template<class _Uty1 = _Ty1,
		class _Uty2 = _Ty2,
		class = enable_if_t<is_copy_constructible<_Uty1>::value
						&& is_copy_constructible<_Uty2>::value>,
		enable_if_t<!is_convertible<const _Uty1&, _Uty1>::value
				|| !is_convertible<const _Uty2&, _Uty2>::value, int> = 0>
		constexpr explicit pair(const _Ty1& _Val1, const _Ty2& _Val2)
		: first(_Val1), second(_Val2)
		{	
		}

	pair(const pair&) = default;
	pair(pair&&) = default;

	template<class _Other1,
		class _Other2,
		class = enable_if_t<is_constructible<_Ty1, const _Other1&>::value
						&& is_constructible<_Ty2, const _Other2&>::value>,
		enable_if_t<is_convertible<const _Other1&, _Ty1>::value
				&& is_convertible<const _Other2&, _Ty2>::value, int> = 0>
		constexpr pair(const pair<_Other1, _Other2>& _Right)
		: first(_Right.first), second(_Right.second)
		{	
		}

	template<class _Other1,
		class _Other2,
		class = enable_if_t<is_constructible<_Ty1, const _Other1&>::value
						&& is_constructible<_Ty2, const _Other2&>::value>,
		enable_if_t<!is_convertible<const _Other1&, _Ty1>::value
				|| !is_convertible<const _Other2&, _Ty2>::value, int> = 0>
		constexpr explicit pair(const pair<_Other1, _Other2>& _Right)
		: first(_Right.first), second(_Right.second)
		{	
		}

	template<class _Other1,
		class _Other2>
		_Myt& operator=(const pair<_Other1, _Other2>& _Right)
		{	
		first = _Right.first;
		second = _Right.second;
		return (*this);
		}

	template<class _Tuple1,
		class _Tuple2,
		size_t... _Indexes1,
		size_t... _Indexes2> inline
		pair(_Tuple1& _Val1,
			_Tuple2& _Val2,
			integer_sequence<size_t, _Indexes1...>,
			integer_sequence<size_t, _Indexes2...>);

	template<class... _Types1,
		class... _Types2> inline
		pair(piecewise_construct_t,
			tuple<_Types1...> _Val1,
			tuple<_Types2...> _Val2);

	template<class _Other1,
		class _Other2,
		class = enable_if_t<is_constructible<_Ty1, _Other1>::value
						&& is_constructible<_Ty2, _Other2>::value>,
		enable_if_t<is_convertible<_Other1, _Ty1>::value
				&& is_convertible<_Other2, _Ty2>::value, int> = 0>
		constexpr pair(_Other1&& _Val1, _Other2&& _Val2)
			noexcept((is_nothrow_constructible<_Ty1, _Other1>::value && is_nothrow_constructible<_Ty2, _Other2>::value))
		: first(::std:: forward<_Other1>(_Val1)),
				second(::std:: forward<_Other2>(_Val2))
		{	
		}

	template<class _Other1,
		class _Other2,
		class = enable_if_t<is_constructible<_Ty1, _Other1>::value
						&& is_constructible<_Ty2, _Other2>::value>,
		enable_if_t<!is_convertible<_Other1, _Ty1>::value
				|| !is_convertible<_Other2, _Ty2>::value, int> = 0>
		constexpr explicit pair(_Other1&& _Val1, _Other2&& _Val2)
			noexcept((is_nothrow_constructible<_Ty1, _Other1>::value && is_nothrow_constructible<_Ty2, _Other2>::value))
		: first(::std:: forward<_Other1>(_Val1)),
				second(::std:: forward<_Other2>(_Val2))
		{	
		}

	template<class _Other1,
		class _Other2,
		class = enable_if_t<is_constructible<_Ty1, _Other1>::value
						&& is_constructible<_Ty2, _Other2>::value>,
		enable_if_t<is_convertible<_Other1, _Ty1>::value
				&& is_convertible<_Other2, _Ty2>::value, int> = 0>
		constexpr pair(pair<_Other1, _Other2>&& _Right)
			noexcept((is_nothrow_constructible<_Ty1, _Other1>::value && is_nothrow_constructible<_Ty2, _Other2>::value))
		: first(::std:: forward<_Other1>(_Right.first)),
			second(::std:: forward<_Other2>(_Right.second))
		{	
		}

	template<class _Other1,
		class _Other2,
		class = enable_if_t<is_constructible<_Ty1, _Other1>::value
						&& is_constructible<_Ty2, _Other2>::value>,
		enable_if_t<!is_convertible<_Other1, _Ty1>::value
				|| !is_convertible<_Other2, _Ty2>::value, int> = 0>
		constexpr explicit pair(pair<_Other1, _Other2>&& _Right)
			noexcept((is_nothrow_constructible<_Ty1, _Other1>::value && is_nothrow_constructible<_Ty2, _Other2>::value))
		: first(::std:: forward<_Other1>(_Right.first)),
			second(::std:: forward<_Other2>(_Right.second))
		{	
		}

	template<class _Other1,
		class _Other2>
		_Myt& operator=(pair<_Other1, _Other2>&& _Right)
		{	
		first = ::std:: forward<_Other1>(_Right.first);
		second = ::std:: forward<_Other2>(_Right.second);
		return (*this);
		}

	_Myt& operator=(_Myt&& _Right)
		noexcept((is_nothrow_move_assignable<_Ty1>::value && is_nothrow_move_assignable<_Ty2>::value))
		{	
		first = ::std:: forward<_Ty1>(_Right.first);
		second = ::std:: forward<_Ty2>(_Right.second);
		return (*this);
		}

	_Myt& operator=(const _Myt& _Right)
		{	
		first = _Right.first;
		second = _Right.second;
		return (*this);
		}

	_Ty1 first;		
	_Ty2 second;	

	void swap(_Myt& _Right)
		noexcept(_Is_nothrow_swappable<_Ty1>::value && _Is_nothrow_swappable<_Ty2>::value)
		{	
		if (this != &_Right)
			{	
			_Swap_adl(first, _Right.first);
			_Swap_adl(second, _Right.second);
			}
		}
	};

		

template<class _Ty1,
	class _Ty2,
	class = enable_if_t<_Is_swappable<_Ty1>::value && _Is_swappable<_Ty2>::value>> inline
	void swap(pair<_Ty1, _Ty2>& _Left, pair<_Ty1, _Ty2>& _Right)
		noexcept(noexcept(_Left.swap(_Right)))
	{	
	_Left.swap(_Right);
	}

template<class _Ty1,
	class _Ty2> inline
	constexpr bool operator==(const pair<_Ty1, _Ty2>& _Left,
		const pair<_Ty1, _Ty2>& _Right)
	{	
	return (_Left.first == _Right.first && _Left.second == _Right.second);
	}

template<class _Ty1,
	class _Ty2> inline
	constexpr bool operator!=(const pair<_Ty1, _Ty2>& _Left,
		const pair<_Ty1, _Ty2>& _Right)
	{	
	return (!(_Left == _Right));
	}

template<class _Ty1,
	class _Ty2> inline
	constexpr bool operator<(const pair<_Ty1, _Ty2>& _Left,
		const pair<_Ty1, _Ty2>& _Right)
	{	
	return (_Left.first < _Right.first ||
		(!(_Right.first < _Left.first) && _Left.second < _Right.second));
	}

template<class _Ty1,
	class _Ty2> inline
	constexpr bool operator>(const pair<_Ty1, _Ty2>& _Left,
		const pair<_Ty1, _Ty2>& _Right)
	{	
	return (_Right < _Left);
	}

template<class _Ty1,
	class _Ty2> inline
	constexpr bool operator<=(const pair<_Ty1, _Ty2>& _Left,
		const pair<_Ty1, _Ty2>& _Right)
	{	
	return (!(_Right < _Left));
	}

template<class _Ty1,
	class _Ty2> inline
	constexpr bool operator>=(const pair<_Ty1, _Ty2>& _Left,
		const pair<_Ty1, _Ty2>& _Right)
	{	
	return (!(_Left < _Right));
	}

	

template<class _Ty1,
	class _Ty2> inline
	constexpr pair<typename _Unrefwrap<_Ty1>::type,
		typename _Unrefwrap<_Ty2>::type>
		make_pair(_Ty1&& _Val1, _Ty2&& _Val2)
	{	
	typedef pair<typename _Unrefwrap<_Ty1>::type,
		typename _Unrefwrap<_Ty2>::type> _Mypair;
	return (_Mypair(::std:: forward<_Ty1>(_Val1),
		::std:: forward<_Ty2>(_Val2)));
	}

		
	namespace rel_ops
		{	
template<class _Ty> inline
	bool operator!=(const _Ty& _Left, const _Ty& _Right)
	{	
	return (!(_Left == _Right));
	}

template<class _Ty> inline
	bool operator>(const _Ty& _Left, const _Ty& _Right)
	{	
	return (_Right < _Left);
	}

template<class _Ty> inline
	bool operator<=(const _Ty& _Left, const _Ty& _Right)
	{	
	return (!(_Right < _Left));
	}

template<class _Ty> inline
	bool operator>=(const _Ty& _Left, const _Ty& _Right)
	{	
	return (!(_Left < _Right));
	}
		}
}

namespace std {
template<class _Ty,
	size_t _Size>
	class array;

	
template<class _Tuple>
	struct tuple_size;

template<class _Ty,
	size_t _Size>
	struct tuple_size<array<_Ty, _Size> >
		: integral_constant<size_t, _Size>
	{	
	};

template<class _Ty1,
	class _Ty2>
	struct tuple_size<pair<_Ty1, _Ty2> >
	: integral_constant<size_t, 2>
	{	
	};

template<class... _Types>
	struct tuple_size<tuple<_Types...> >
	: integral_constant<size_t, sizeof...(_Types)>
	{	
	};


template<class _Tuple>
	struct tuple_size<const _Tuple>
	: tuple_size<_Tuple>
	{	
	};

template<class _Tuple>
	struct tuple_size<volatile _Tuple>
	: tuple_size<_Tuple>
	{	
	};

template<class _Tuple>
	struct tuple_size<const volatile _Tuple>
	: tuple_size<_Tuple>
	{	
	};

 
template<class _Ty>
	constexpr size_t tuple_size_v = tuple_size<_Ty>::value;
 

	
template<size_t _Index,
	class _Tuple>
	struct tuple_element;

template<size_t _Idx,
	class _Ty,
	size_t _Size>
	struct tuple_element<_Idx, array<_Ty, _Size> >
	{	
	static_assert(_Idx < _Size, "array index out of bounds");

	typedef _Ty type;
	};

template<class _Ty1,
	class _Ty2>
	struct tuple_element<0, pair<_Ty1, _Ty2> >
	{	
	typedef _Ty1 type;
	};

template<class _Ty1,
	class _Ty2>
	struct tuple_element<1, pair<_Ty1, _Ty2> >
	{	
	typedef _Ty2 type;
	};

template<size_t _Index>
	struct tuple_element<_Index, tuple<> >
	{	
	static_assert(_Always_false<integral_constant<size_t, _Index> >::value,
		"tuple index out of bounds");
	};

template<class _This,
	class... _Rest>
	struct tuple_element<0, tuple<_This, _Rest...> >
	{	
	typedef _This type;
	typedef tuple<_This, _Rest...> _Ttype;
	};

template<size_t _Index,
	class _This,
	class... _Rest>
	struct tuple_element<_Index, tuple<_This, _Rest...> >
		: public tuple_element<_Index - 1, tuple<_Rest...> >
	{	
	};


template<size_t _Index,
	class _Tuple>
	struct tuple_element<_Index, const _Tuple>
		: public tuple_element<_Index, _Tuple>
	{	
	typedef tuple_element<_Index, _Tuple> _Mybase;
	typedef typename add_const<typename _Mybase::type>::type type;
	};

template<size_t _Index,
	class _Tuple>
	struct tuple_element<_Index, volatile _Tuple>
		: public tuple_element<_Index, _Tuple>
	{	
	typedef tuple_element<_Index, _Tuple> _Mybase;
	typedef typename add_volatile<typename _Mybase::type>::type type;
	};

template<size_t _Index,
	class _Tuple>
	struct tuple_element<_Index, const volatile _Tuple>
		: public tuple_element<_Index, _Tuple>
	{	
	typedef tuple_element<_Index, _Tuple> _Mybase;
	typedef typename add_cv<typename _Mybase::type>::type type;
	};

template<size_t _Index,
	class _Tuple>
	using tuple_element_t = typename tuple_element<_Index, _Tuple>::type;

	
template<class _Ret,
	class _Pair> inline
	constexpr _Ret _Pair_get(_Pair& _Pr,
		integral_constant<size_t, 0>) noexcept
	{	
	return (_Pr.first);
	}

template<class _Ret,
	class _Pair> inline
	constexpr _Ret _Pair_get(_Pair& _Pr,
		integral_constant<size_t, 1>) noexcept
	{	
	return (_Pr.second);
	}

template<size_t _Idx,
	class _Ty1,
	class _Ty2> inline
	constexpr typename tuple_element<_Idx, pair<_Ty1, _Ty2> >::type&
		get(pair<_Ty1, _Ty2>& _Pr) noexcept
	{	
	typedef typename tuple_element<_Idx, pair<_Ty1, _Ty2> >::type& _Rtype;
	return (_Pair_get<_Rtype>(_Pr, integral_constant<size_t, _Idx>()));
	}

template<class _Ty1,
	class _Ty2> inline
	constexpr _Ty1& get(pair<_Ty1, _Ty2>& _Pr) noexcept
	{	
	return (::std:: get<0>(_Pr));
	}

template<class _Ty2,
	class _Ty1> inline
	constexpr _Ty2& get(pair<_Ty1, _Ty2>& _Pr) noexcept
	{	
	return (::std:: get<1>(_Pr));
	}

template<size_t _Idx,
	class _Ty1,
	class _Ty2> inline
	constexpr const typename tuple_element<_Idx, pair<_Ty1, _Ty2> >::type&
		get(const pair<_Ty1, _Ty2>& _Pr) noexcept
	{	
	typedef const typename tuple_element<_Idx, pair<_Ty1, _Ty2> >::type&
		_Ctype;
	return (_Pair_get<_Ctype>(_Pr, integral_constant<size_t, _Idx>()));
	}

template<class _Ty1,
	class _Ty2> inline
	constexpr const _Ty1& get(const pair<_Ty1, _Ty2>& _Pr) noexcept
	{	
	return (::std:: get<0>(_Pr));
	}

template<class _Ty2,
	class _Ty1> inline
	constexpr const _Ty2& get(const pair<_Ty1, _Ty2>& _Pr) noexcept
	{	
	return (::std:: get<1>(_Pr));
	}

template<size_t _Idx,
	class _Ty1,
	class _Ty2> inline
	constexpr typename tuple_element<_Idx, pair<_Ty1, _Ty2> >::type&&
		get(pair<_Ty1, _Ty2>&& _Pr) noexcept
	{	
	typedef typename tuple_element<_Idx, pair<_Ty1, _Ty2> >::type&& _RRtype;
	return (::std:: forward<_RRtype>(::std:: get<_Idx>(_Pr)));
	}

template<class _Ty1,
	class _Ty2> inline
	constexpr _Ty1&& get(pair<_Ty1, _Ty2>&& _Pr) noexcept
	{	
	return (::std:: get<0>(::std:: move(_Pr)));
	}

template<class _Ty2,
	class _Ty1> inline
	constexpr _Ty2&& get(pair<_Ty1, _Ty2>&& _Pr) noexcept
	{	
	return (::std:: get<1>(::std:: move(_Pr)));
	}

	
template<class _Ty,
	class _Other = _Ty> inline
	_Ty exchange(_Ty& _Val, _Other&& _New_val)
	{	
	_Ty _Old_val = ::std:: move(_Val);
	_Val = ::std:: forward<_Other>(_New_val);
	return (_Old_val);
	}

	
template<class _Ty> inline
	constexpr add_const_t<_Ty>& as_const(_Ty& _Val) noexcept
	{	
	return (_Val);
	}

template<class _Ty>
	void as_const(const _Ty&&) = delete;
}


namespace std {
namespace tr1 {	
using ::std:: get;
using ::std:: tuple_element;
using ::std:: tuple_size;
}	
}


 
 #pragma warning(pop)
 #pragma pack(pop)









 #pragma pack(push,8)
 #pragma warning(push,3)
 
 

namespace std {
		

 













  
 

		
struct _Container_proxy;
struct _Container_base12;
struct _Iterator_base12;

struct _Container_base0
	{	
	void _Orphan_all()
		{	
		}

	void _Swap_all(_Container_base0&)
		{	
		}
	};

struct _Iterator_base0
	{	
	void _Adopt(const void *)
		{	
		}

	const _Container_base0 *_Getcont() const
		{	
		return (0);
		}
	};

		
struct _Container_proxy
	{	
	_Container_proxy()
		: _Mycont(0), _Myfirstiter(0)
		{	
		}

	const _Container_base12 *_Mycont;
	_Iterator_base12 *_Myfirstiter;
	};

struct _Container_base12
	{	
public:
	_Container_base12()
		: _Myproxy(0)
		{	
		}

	_Container_base12(const _Container_base12&)
		: _Myproxy(0)
		{	
		}

	_Container_base12& operator=(const _Container_base12&)
		{	
		return (*this);
		}

	~_Container_base12() noexcept
		{	
		_Orphan_all();
		}

	_Iterator_base12 **_Getpfirst() const
		{	
		return (_Myproxy == 0 ? 0 : &_Myproxy->_Myfirstiter);
		}

	void _Orphan_all();	
	void _Swap_all(_Container_base12&);	

	_Container_proxy *_Myproxy;
	};

struct _Iterator_base12
	{	
public:
	_Iterator_base12()
		: _Myproxy(0), _Mynextiter(0)
		{	
		}

	_Iterator_base12(const _Iterator_base12& _Right)
		: _Myproxy(0), _Mynextiter(0)
		{	
		*this = _Right;
		}

	_Iterator_base12& operator=(const _Iterator_base12& _Right)
		{	
		if (_Myproxy == _Right._Myproxy)
			;
		else if (_Right._Myproxy != 0)
			_Adopt(_Right._Myproxy->_Mycont);
		else
			{	
 



			}

		return (*this);
		}

	~_Iterator_base12() noexcept
		{	
 



		}

	void _Adopt(const _Container_base12 *_Parent)
		{	
		if (_Parent == 0)
			{	
 



			}
		else
			{	
			_Container_proxy *_Parent_proxy = _Parent->_Myproxy;

 










			_Myproxy = _Parent_proxy;
 
			}
		}

	void _Clrcont()
		{	
		_Myproxy = 0;
		}

	const _Container_base12 *_Getcont() const
		{	
		return (_Myproxy == 0 ? 0 : _Myproxy->_Mycont);
		}

	_Iterator_base12 **_Getpnext()
		{	
		return (&_Mynextiter);
		}

	void _Orphan_me()
		{	
 












		}

	_Container_proxy *_Myproxy;
	_Iterator_base12 *_Mynextiter;
	};

		
inline void _Container_base12::_Orphan_all()
	{	
 










	}

inline void _Container_base12::_Swap_all(_Container_base12& _Right)
	{	
 



	_Container_proxy *_Temp = _Myproxy;
	_Myproxy = _Right._Myproxy;
	_Right._Myproxy = _Temp;

	if (_Myproxy != 0)
		_Myproxy->_Mycont = (_Container_base12 *)this;
	if (_Right._Myproxy != 0)
		_Right._Myproxy->_Mycont = (_Container_base12 *)&_Right;
	}

 
typedef _Container_base0 _Container_base;
typedef _Iterator_base0 _Iterator_base;

 




	
struct _Zero_then_variadic_args_t
	{	
	};	

struct _One_then_variadic_args_t
	{	
	};	

template<class _Ty1,
	class _Ty2,
	bool = is_empty<_Ty1>::value && !is_final<_Ty1>::value>
	class _Compressed_pair final
		: private _Ty1

	{	
private:
	_Ty2 _Myval2;

	typedef _Ty1 _Mybase;	

public:
	template<class... _Other2>
		constexpr explicit _Compressed_pair(_Zero_then_variadic_args_t,
			_Other2&&... _Val2)
		: _Ty1(), _Myval2(::std:: forward<_Other2>(_Val2)...)
		{	
		}

	template<class _Other1,
		class... _Other2>
		_Compressed_pair(_One_then_variadic_args_t,
			_Other1&& _Val1, _Other2&&... _Val2)
		: _Ty1(::std:: forward<_Other1>(_Val1)),
			_Myval2(::std:: forward<_Other2>(_Val2)...)
		{	
		}


	_Ty1& _Get_first() noexcept
		{	
		return (*this);
		}

	const _Ty1& _Get_first() const noexcept
		{	
		return (*this);
		}

	volatile _Ty1& _Get_first() volatile noexcept
		{	
		return (*this);
		}

	const volatile _Ty1& _Get_first() const volatile noexcept
		{	
		return (*this);
		}

	_Ty2& _Get_second() noexcept
		{	
		return (_Myval2);
		}

	const _Ty2& _Get_second() const noexcept
		{	
		return (_Myval2);
		}

	volatile _Ty2& _Get_second() volatile noexcept
		{	
		return (_Myval2);
		}

	const volatile _Ty2& _Get_second() const volatile noexcept
		{	
		return (_Myval2);
		}
	};

template<class _Ty1,
	class _Ty2>
	class _Compressed_pair<_Ty1, _Ty2, false> final

	{	
private:
	_Ty1 _Myval1;
	_Ty2 _Myval2;

public:
	template<class... _Other2>
		constexpr explicit _Compressed_pair(_Zero_then_variadic_args_t,
			_Other2&&... _Val2)
		: _Myval1(), _Myval2(::std:: forward<_Other2>(_Val2)...)
		{	
		}

	template<class _Other1,
		class... _Other2>
		_Compressed_pair(_One_then_variadic_args_t,
			_Other1&& _Val1, _Other2&&... _Val2)
		: _Myval1(::std:: forward<_Other1>(_Val1)),
			_Myval2(::std:: forward<_Other2>(_Val2)...)
		{	
		}


	_Ty1& _Get_first() noexcept
		{	
		return (_Myval1);
		}

	const _Ty1& _Get_first() const noexcept
		{	
		return (_Myval1);
		}

	volatile _Ty1& _Get_first() volatile noexcept
		{	
		return (_Myval1);
		}

	const volatile _Ty1& _Get_first() const volatile noexcept
		{	
		return (_Myval1);
		}

	_Ty2& _Get_second() noexcept
		{	
		return (_Myval2);
		}

	const _Ty2& _Get_second() const noexcept
		{	
		return (_Myval2);
		}

	volatile _Ty2& _Get_second() volatile noexcept
		{	
		return (_Myval2);
		}

	const volatile _Ty2& _Get_second() const volatile noexcept
		{	
		return (_Myval2);
		}
	};

		
template<class _Ty,
	class = void>
	struct _Is_checked_helper
		: false_type
	{	
	};

template<class _Ty>
	struct _Is_checked_helper<_Ty, void_t<
		typename _Ty::_Unchecked_type> >
		: true_type
	{	
	};

		
template<class _Iter> inline
	typename _Is_checked_helper<_Iter>::type _Is_checked(_Iter)
	{	
	return {};
	}

		
template<class _Iter> inline
	_Iter _Unchecked(_Iter _Src)
	{	
	return (_Src);
	}

 


		
 

template<class _Iter> inline
	decltype(_Unchecked(::std:: declval<_Iter>())) _Unchecked_idl0(_Iter _Src)
	{	
	return (_Unchecked(_Src));
	}

 









		
template<class _Iter,
	class _UIter> inline
	_Iter& _Rechecked(_Iter& _Dest, _UIter _Src)
	{	
	_Dest = _Src;
	return (_Dest);
	}

		




















		
		
		


















		
		
struct input_iterator_tag
	{	
	};

struct _Mutable_iterator_tag	
	{	
	};

struct output_iterator_tag
	: _Mutable_iterator_tag
	{	
	};

struct forward_iterator_tag
	: input_iterator_tag, _Mutable_iterator_tag
	{	
	};

struct bidirectional_iterator_tag
	: forward_iterator_tag
	{	
	};

struct random_access_iterator_tag
	: bidirectional_iterator_tag
	{	
	};

		
struct _General_ptr_iterator_tag
	{	
	};

struct _Trivially_copyable_ptr_iterator_tag
	: _General_ptr_iterator_tag
	{	
	};

struct _Really_trivial_ptr_iterator_tag
	: _Trivially_copyable_ptr_iterator_tag
	{	
	};

	
struct _Any_tag
	{	
	constexpr _Any_tag() noexcept = default;
	template<class _Ty>
		constexpr _Any_tag(_Ty&&) noexcept {}
	};

		
template<class _Category,
	class _Ty,
	class _Diff = ptrdiff_t,
	class _Pointer = _Ty *,
	class _Reference = _Ty&>
	struct iterator
	{	
	typedef _Category iterator_category;
	typedef _Ty value_type;
	typedef _Diff difference_type;

	typedef _Pointer pointer;
	typedef _Reference reference;
	};

template<class _Category,
	class _Ty,
	class _Diff,
	class _Pointer,
	class _Reference,
	class _Base>
	struct _Iterator012
		: public _Base
	{	
	typedef _Category iterator_category;
	typedef _Ty value_type;
	typedef _Diff difference_type;

	typedef _Pointer pointer;
	typedef _Reference reference;
	};


typedef iterator<output_iterator_tag, void, void, void, void> _Outit;

		
template<class,
	class = void>
	struct _Iterator_traits_base
	{	
	};

template<class _Iter>
	struct _Iterator_traits_base<_Iter, void_t<
		typename _Iter::iterator_category,
		typename _Iter::value_type,
		typename _Iter::difference_type,
		typename _Iter::pointer,
		typename _Iter::reference
		> >
	{	
	typedef typename _Iter::iterator_category iterator_category;
	typedef typename _Iter::value_type value_type;
	typedef typename _Iter::difference_type difference_type;

	typedef typename _Iter::pointer pointer;
	typedef typename _Iter::reference reference;
	};

template<class _Iter>
	struct iterator_traits
		: _Iterator_traits_base<_Iter>
	{	
	};

template<class _Ty>
	struct iterator_traits<_Ty *>
	{	
	typedef random_access_iterator_tag iterator_category;
	typedef _Ty value_type;
	typedef ptrdiff_t difference_type;

	typedef _Ty *pointer;
	typedef _Ty& reference;
	};

template<class _Ty>
	struct iterator_traits<const _Ty *>
	{	
	typedef random_access_iterator_tag iterator_category;
	typedef _Ty value_type;
	typedef ptrdiff_t difference_type;

	typedef const _Ty *pointer;
	typedef const _Ty& reference;
	};

		
template<class _Iter>
	using _Iter_value_t = typename iterator_traits<_Iter>::value_type;

		
template<class _Iter>
	using _Iter_diff_t = typename iterator_traits<_Iter>::difference_type;

		
template<class _Iter>
	using _Iter_cat_t = typename iterator_traits<_Iter>::iterator_category;

		
template<class _Ty,
	class = void>
	struct _Is_iterator
		: false_type
	{	
	};

template<class _Ty>
	struct _Is_iterator<_Ty, void_t<
		typename iterator_traits<_Ty>::iterator_category
		> >
		: true_type
	{	
	};


		
 
template<class _Iter,
	class _Diff> inline
	auto _Unchecked_n(_Iter _Src, _Diff)
	{	
	return (_Unchecked(_Src));
	}
 




























		
template<class _Ty1,
	class _Ty2>
	struct _Is_same_size
		: bool_constant<sizeof(_Ty1) == sizeof(_Ty2)>
	{	
	};

		
template<class _Elem,
	bool _Is_enum = is_enum<_Elem>::value>
	struct _Unwrap_enum
	{	
	typedef underlying_type_t<_Elem> type;
	};

template<class _Elem>
	struct _Unwrap_enum<_Elem, false>
	{	
	typedef _Elem type;
	};

		
template<class _Ty1,
	class _Ty2>
	struct _Both_or_neither_bool
		: bool_constant<is_same<bool, _Ty1>::value == is_same<bool, _Ty2>::value>
	{	
	};

		
template<class _Source,
	class _Dest>
	struct _Ptr_cat_helper
	{	
	typedef typename _Unwrap_enum<_Source>::type _USource;
	typedef typename _Unwrap_enum<_Dest>::type _UDest;
	typedef conditional_t<
		conjunction<
			_Is_same_size<_USource, _UDest>,
			is_integral<_USource>,
			is_integral<_UDest>,
			_Both_or_neither_bool<_USource, _UDest>,
			
			negation<is_volatile<_Source>>,
			negation<is_volatile<_Dest>>
		>::value,
		_Really_trivial_ptr_iterator_tag,
		_General_ptr_iterator_tag> type;
	};

template<class _Elem>
	struct _Ptr_cat_helper<_Elem, _Elem>
	{	
	typedef conditional_t<
		is_trivially_copyable<_Elem>::value,
		conditional_t<is_trivial<_Elem>::value,
			_Really_trivial_ptr_iterator_tag,
			_Trivially_copyable_ptr_iterator_tag>,
		_General_ptr_iterator_tag> type;
	};

template<class _Anything>
	struct _Ptr_cat_helper<_Anything *, const _Anything *>
	{	
	typedef _Really_trivial_ptr_iterator_tag type;
	};

template<class _Source,
	class _Dest> inline
	_General_ptr_iterator_tag _Ptr_copy_cat(const _Source&, const _Dest&)
	{	
	return {};
	}

template<class _Source,
	class _Dest> inline
	conditional_t<is_trivially_assignable<_Dest&, _Source&>::value,
		typename _Ptr_cat_helper<remove_const_t<_Source>, _Dest>::type,
		_General_ptr_iterator_tag>
		_Ptr_copy_cat(_Source * const&, _Dest * const&)
	{	
	return {};
	}

template<class _Source,
	class _Dest> inline
	_General_ptr_iterator_tag _Ptr_move_cat(const _Source&, const _Dest&)
	{	
	return {};
	}

template<class _Source,
	class _Dest> inline
	conditional_t<is_trivially_assignable<_Dest&, _Source>::value,
		typename _Ptr_cat_helper<remove_const_t<_Source>, _Dest>::type,
		_General_ptr_iterator_tag>
		_Ptr_move_cat(_Source * const&, _Dest * const&)
	{	
	return {};
	}

		

 
  
  
  
  
  
  
  
  
  
  

 











































































































































































































 
  
 





















		
		
template<class _InIt,
	class _Diff> inline
	void _Advance1(_InIt& _Where, _Diff _Off, input_iterator_tag)
	{	
 




	for (; 0 < _Off; --_Off)
		++_Where;
	}

template<class _BidIt,
	class _Diff> inline
	void _Advance1(_BidIt& _Where, _Diff _Off, bidirectional_iterator_tag)
	{	
	for (; 0 < _Off; --_Off)
		++_Where;
	for (; _Off < 0; ++_Off)
		--_Where;
	}

template<class _RanIt,
	class _Diff> inline
	void _Advance1(_RanIt& _Where, _Diff _Off, random_access_iterator_tag)
	{	
	_Where += _Off;
	}

template<class _InIt,
	class _Diff> inline
	void advance(_InIt& _Where, _Diff _Off)
	{	
		
	_Advance1(_Where, _Off, _Iter_cat_t<remove_const_t<_InIt>>());
	}

		
template<class _InIt> inline
	_Iter_diff_t<_InIt>
		_Distance1(_InIt _First, _InIt _Last, input_iterator_tag)
	{	
	_Iter_diff_t<_InIt> _Off = 0;
	for (; _First != _Last; ++_First)
		++_Off;

	return (_Off);
	}

template<class _RanIt> inline
	_Iter_diff_t<_RanIt>
		_Distance1(_RanIt _First, _RanIt _Last, random_access_iterator_tag)
	{	
 







	return (_Last - _First);
	}

template<class _InIt> inline
	_Iter_diff_t<_InIt>
		distance(_InIt _First, _InIt _Last)
	{	
	return (_Distance1(_First, _Last, _Iter_cat_t<_InIt>()));
	}

		
template<class _InIt> inline
	_InIt next(_InIt _First, _Iter_diff_t<_InIt> _Off = 1)
	{	
	static_assert(is_base_of<input_iterator_tag,
		typename iterator_traits<_InIt>::iterator_category>::value,
		"next requires input iterator");

	::std:: advance(_First, _Off);
	return (_First);
	}

		
template<class _BidIt> inline
	_BidIt prev(_BidIt _First, _Iter_diff_t<_BidIt> _Off = 1)
	{	
	static_assert(is_base_of<bidirectional_iterator_tag,
		typename iterator_traits<_BidIt>::iterator_category>::value,
		"prev requires bidirectional iterator");

	::std:: advance(_First, -_Off);
	return (_First);
	}

		
template<class _Ty>
	struct pointer_traits;

template<class _RanIt>
	class reverse_iterator
		: public iterator<
			typename iterator_traits<_RanIt>::iterator_category,
			typename iterator_traits<_RanIt>::value_type,
			typename iterator_traits<_RanIt>::difference_type,
			typename iterator_traits<_RanIt>::pointer,
			typename iterator_traits<_RanIt>::reference>
	{	
	typedef reverse_iterator<_RanIt> _Myt;

public:
	typedef typename iterator_traits<_RanIt>::difference_type difference_type;
	typedef typename iterator_traits<_RanIt>::pointer pointer;
	typedef typename iterator_traits<_RanIt>::reference reference;
	typedef _RanIt iterator_type;

	reverse_iterator()
		: current()
		{	
		}

	explicit reverse_iterator(_RanIt _Right)
		: current(_Right)
		{	
		}

	template<class _Other>
		reverse_iterator(const reverse_iterator<_Other>& _Right)
		: current(_Right.base())
		{	
		}

	template<class _Other>
		_Myt& operator=(const reverse_iterator<_Other>& _Right)
		{	
		current = _Right.base();
		return (*this);
		}

	_RanIt base() const
		{	
		return (current);
		}

	reference operator*() const
		{	
		_RanIt _Tmp = current;
		return (*--_Tmp);
		}

	pointer operator->() const
		{	
		return (pointer_traits<pointer>::pointer_to(**this));
		}

	_Myt& operator++()
		{	
		--current;
		return (*this);
		}

	_Myt operator++(int)
		{	
		_Myt _Tmp = *this;
		--current;
		return (_Tmp);
		}

	_Myt& operator--()
		{	
		++current;
		return (*this);
		}

	_Myt operator--(int)
		{	
		_Myt _Tmp = *this;
		++current;
		return (_Tmp);
		}



	_Myt& operator+=(difference_type _Off)
		{	
		current -= _Off;
		return (*this);
		}

	_Myt operator+(difference_type _Off) const
		{	
		return (_Myt(current - _Off));
		}

	_Myt& operator-=(difference_type _Off)
		{	
		current += _Off;
		return (*this);
		}

	_Myt operator-(difference_type _Off) const
		{	
		return (_Myt(current + _Off));
		}

	reference operator[](difference_type _Off) const
		{	
		return (*(*this + _Off));
		}

protected:
	_RanIt current;	
	};

template<class _RanIt>
	struct _Is_checked_helper<reverse_iterator<_RanIt> >
		: public _Is_checked_helper<_RanIt>
	{	
	};

		
template<class _RanIt> inline
	reverse_iterator<_RanIt> operator+(
		typename reverse_iterator<_RanIt>::difference_type _Off,
		const reverse_iterator<_RanIt>& _Right)
	{	
	return (_Right + _Off);
	}

template<class _RanIt1,
	class _RanIt2>
	auto inline operator-(const reverse_iterator<_RanIt1>& _Left,
		const reverse_iterator<_RanIt2>& _Right)
			-> decltype(_Right.base() - _Left.base())
	{	
	return (_Right.base() - _Left.base());
	}

template<class _RanIt1,
	class _RanIt2> inline
	bool operator==(const reverse_iterator<_RanIt1>& _Left,
		const reverse_iterator<_RanIt2>& _Right)
	{	
	return (_Left.base() == _Right.base());
	}

template<class _RanIt1,
	class _RanIt2> inline
	bool operator!=(const reverse_iterator<_RanIt1>& _Left,
		const reverse_iterator<_RanIt2>& _Right)
	{	
	return (!(_Left == _Right));
	}

template<class _RanIt1,
	class _RanIt2> inline
	bool operator<(const reverse_iterator<_RanIt1>& _Left,
		const reverse_iterator<_RanIt2>& _Right)
	{	
	return (_Right.base() < _Left.base());
	}

template<class _RanIt1,
	class _RanIt2> inline
	bool operator>(const reverse_iterator<_RanIt1>& _Left,
		const reverse_iterator<_RanIt2>& _Right)
	{	
	return (_Right < _Left);
	}

template<class _RanIt1,
	class _RanIt2> inline
	bool operator<=(const reverse_iterator<_RanIt1>& _Left,
		const reverse_iterator<_RanIt2>& _Right)
	{	
	return (!(_Right < _Left));
	}

template<class _RanIt1,
	class _RanIt2> inline
	bool operator>=(const reverse_iterator<_RanIt1>& _Left,
		const reverse_iterator<_RanIt2>& _Right)
	{	
	return (!(_Left < _Right));
	}

		
template<class _RanIt> inline
	reverse_iterator<_RanIt> make_reverse_iterator(_RanIt _Iter)
	{	
	return (reverse_iterator<_RanIt>(_Iter));
	}

		

template<class _Container>
	auto inline begin(_Container& _Cont) -> decltype(_Cont.begin())
	{	
	return (_Cont.begin());
	}

template<class _Container>
	auto inline begin(const _Container& _Cont) -> decltype(_Cont.begin())
	{	
	return (_Cont.begin());
	}

template<class _Container>
	auto inline end(_Container& _Cont) -> decltype(_Cont.end())
	{	
	return (_Cont.end());
	}

template<class _Container>
	auto inline end(const _Container& _Cont) -> decltype(_Cont.end())
	{	
	return (_Cont.end());
	}

template<class _Ty,
	size_t _Size> inline
	constexpr _Ty *begin(_Ty (&_Array)[_Size]) noexcept
	{	
	return (_Array);
	}

template<class _Ty,
	size_t _Size> inline
	constexpr _Ty *end(_Ty (&_Array)[_Size]) noexcept
	{	
	return (_Array + _Size);
	}

		
template<class _Container>
	constexpr auto inline cbegin(const _Container& _Cont)
		noexcept(noexcept(::std:: begin(_Cont)))
		-> decltype(::std:: begin(_Cont))
	{	
	return (::std:: begin(_Cont));
	}

template<class _Container>
	constexpr auto inline cend(const _Container& _Cont)
		noexcept(noexcept(::std:: end(_Cont)))
		-> decltype(::std:: end(_Cont))
	{	
	return (::std:: end(_Cont));
	}

		
template<class _Container>
	auto inline rbegin(_Container& _Cont) -> decltype(_Cont.rbegin())
	{	
	return (_Cont.rbegin());
	}

template<class _Container>
	auto inline rbegin(const _Container& _Cont) -> decltype(_Cont.rbegin())
	{	
	return (_Cont.rbegin());
	}

template<class _Container>
	auto inline rend(_Container& _Cont) -> decltype(_Cont.rend())
	{	
	return (_Cont.rend());
	}

template<class _Container>
	auto inline rend(const _Container& _Cont) -> decltype(_Cont.rend())
	{	
	return (_Cont.rend());
	}

template<class _Ty,
	size_t _Size> inline
	reverse_iterator<_Ty *> rbegin(_Ty (&_Array)[_Size])
	{	
	return (reverse_iterator<_Ty *>(_Array + _Size));
	}

template<class _Ty,
	size_t _Size> inline
	reverse_iterator<_Ty *> rend(_Ty (&_Array)[_Size])
	{	
	return (reverse_iterator<_Ty *>(_Array));
	}

template<class _Elem> inline
	reverse_iterator<const _Elem *>
		rbegin(::std:: initializer_list<_Elem> _Ilist)
	{	
	return (reverse_iterator<const _Elem *>(_Ilist.end()));
	}

template<class _Elem> inline
	reverse_iterator<const _Elem *>
		rend(::std:: initializer_list<_Elem> _Ilist)
	{	
	return (reverse_iterator<const _Elem *>(_Ilist.begin()));
	}

		
template<class _Container>
	auto inline crbegin(const _Container& _Cont)
		-> decltype(::std:: rbegin(_Cont))
	{	
	return (::std:: rbegin(_Cont));
	}

template<class _Container>
	auto inline crend(const _Container& _Cont)
		-> decltype(::std:: rend(_Cont))
	{	
	return (::std:: rend(_Cont));
	}


template<class _Container>
	constexpr auto inline size(const _Container& _Cont)
		-> decltype(_Cont.size())
	{	
	return (_Cont.size());
	}

template<class _Ty,
	size_t _Size> inline
	constexpr size_t size(const _Ty(&)[_Size]) noexcept
	{	
	return (_Size);
	}

template<class _Container>
	constexpr auto inline empty(const _Container& _Cont)
		-> decltype(_Cont.empty())
	{	
	return (_Cont.empty());
	}

template<class _Ty,
	size_t _Size> inline
	constexpr bool empty(const _Ty(&)[_Size]) noexcept
	{	
	return (false);
	}

template<class _Elem> inline
	constexpr bool empty(
		::std:: initializer_list<_Elem> _Ilist) noexcept
	{	
	return (_Ilist.size() == 0);
	}

template<class _Container>
	constexpr auto inline data(_Container& _Cont)
		-> decltype(_Cont.data())
	{	
	return (_Cont.data());
	}

template<class _Container>
	constexpr auto inline data(const _Container& _Cont)
		-> decltype(_Cont.data())
	{	
	return (_Cont.data());
	}

template<class _Ty,
	size_t _Size> inline
	constexpr _Ty *data(_Ty(&_Array)[_Size]) noexcept
	{	
	return (_Array);
	}

template<class _Elem> inline
	constexpr const _Elem *data(
		::std:: initializer_list<_Elem> _Ilist) noexcept
	{	
	return (_Ilist.begin());
	}

		
template<class _Ty,
	size_t _Size>
	class _Array_const_iterator
		: public _Iterator012<random_access_iterator_tag,
			_Ty,
			ptrdiff_t,
			const _Ty *,
			const _Ty&,
			_Iterator_base>
	{	
public:
	typedef _Array_const_iterator<_Ty, _Size> _Myiter;
	typedef random_access_iterator_tag iterator_category;

	typedef _Ty value_type;
	typedef size_t size_type;
	typedef ptrdiff_t difference_type;
	typedef const _Ty *pointer;
	typedef const _Ty& reference;
	enum {_EEN_SIZE = _Size};	
 
	_Array_const_iterator()
		: _Ptr(0)
		{	
		}

	explicit _Array_const_iterator(pointer _Parg, size_t _Off = 0)
		: _Ptr(_Parg + _Off)
		{	
		}

	typedef pointer _Unchecked_type;

	_Myiter& _Rechecked(_Unchecked_type _Right)
		{	
		_Ptr = _Right;
		return (*this);
		}

	_Unchecked_type _Unchecked() const
		{	
		return (_Ptr);
		}

	reference operator*() const
		{	
		return (*_Ptr);
		}

	pointer operator->() const
		{	
		return (pointer_traits<pointer>::pointer_to(**this));
		}

	_Myiter& operator++()
		{	
		++_Ptr;
		return (*this);
		}

	_Myiter operator++(int)
		{	
		_Myiter _Tmp = *this;
		++*this;
		return (_Tmp);
		}

	_Myiter& operator--()
		{	
		--_Ptr;
		return (*this);
		}

	_Myiter operator--(int)
		{	
		_Myiter _Tmp = *this;
		--*this;
		return (_Tmp);
		}

	_Myiter& operator+=(difference_type _Off)
		{	
		_Ptr += _Off;
		return (*this);
		}

	_Myiter operator+(difference_type _Off) const
		{	
		_Myiter _Tmp = *this;
		return (_Tmp += _Off);
		}

	_Myiter& operator-=(difference_type _Off)
		{	
		return (*this += -_Off);
		}

	_Myiter operator-(difference_type _Off) const
		{	
		_Myiter _Tmp = *this;
		return (_Tmp -= _Off);
		}

	difference_type operator-(const _Myiter& _Right) const
		{	
		return (_Ptr - _Right._Ptr);
		}

	reference operator[](difference_type _Off) const
		{	
		return (*(*this + _Off));
		}

	bool operator==(const _Myiter& _Right) const
		{	
		return (_Ptr == _Right._Ptr);
		}

	bool operator!=(const _Myiter& _Right) const
		{	
		return (!(*this == _Right));
		}

	bool operator<(const _Myiter& _Right) const
		{	
		return (_Ptr < _Right._Ptr);
		}

	bool operator>(const _Myiter& _Right) const
		{	
		return (_Right < *this);
		}

	bool operator<=(const _Myiter& _Right) const
		{	
		return (!(_Right < *this));
		}

	bool operator>=(const _Myiter& _Right) const
		{	
		return (!(*this < _Right));
		}

	pointer _Ptr;	

 









































































































































































































	};

template<class _Ty,
	size_t _Size> inline
	typename _Array_const_iterator<_Ty, _Size>::_Unchecked_type
		_Unchecked(_Array_const_iterator<_Ty, _Size> _Iter)
	{	
	return (_Iter._Unchecked());
	}

template<class _Ty,
	size_t _Size> inline
	_Array_const_iterator<_Ty, _Size>&
		_Rechecked(_Array_const_iterator<_Ty, _Size>& _Iter,
			typename _Array_const_iterator<_Ty, _Size>
				::_Unchecked_type _Right)
	{	
	return (_Iter._Rechecked(_Right));
	}

template<class _Ty,
	size_t _Size> inline
	_Array_const_iterator<_Ty, _Size> operator+(
		typename _Array_const_iterator<_Ty, _Size>::difference_type _Off,
		_Array_const_iterator<_Ty, _Size> _Next)
	{	
	return (_Next += _Off);
	}

		
template<class _Ty,
	size_t _Size>
	class _Array_iterator
		: public _Array_const_iterator<_Ty, _Size>
	{	
public:
	typedef _Array_iterator<_Ty, _Size> _Myiter;
	typedef _Array_const_iterator<_Ty, _Size> _Mybase;
	typedef random_access_iterator_tag iterator_category;

	typedef _Ty value_type;
	typedef size_t size_type;
	typedef ptrdiff_t difference_type;
	typedef _Ty *pointer;
	typedef _Ty& reference;

	_Array_iterator()
		{	
		}

	explicit _Array_iterator(pointer _Parg, size_t _Off = 0)
		: _Mybase(_Parg, _Off)
		{	
		}
	enum {_EEN_SIZE = _Size};	
	typedef pointer _Unchecked_type;

	_Myiter& _Rechecked(_Unchecked_type _Right)
		{	
		((_Mybase *)this)->_Rechecked(_Right);
		return (*this);
		}

	_Unchecked_type _Unchecked() const
		{	
		return ((pointer)((_Mybase *)this)->_Unchecked());
		}

	reference operator*() const
		{	
		return ((reference)**(_Mybase *)this);
		}

	pointer operator->() const
		{	
		return (pointer_traits<pointer>::pointer_to(**this));
		}

	_Myiter& operator++()
		{	
		++*(_Mybase *)this;
		return (*this);
		}

	_Myiter operator++(int)
		{	
		_Myiter _Tmp = *this;
		++*this;
		return (_Tmp);
		}

	_Myiter& operator--()
		{	
		--*(_Mybase *)this;
		return (*this);
		}

	_Myiter operator--(int)
		{	
		_Myiter _Tmp = *this;
		--*this;
		return (_Tmp);
		}

	_Myiter& operator+=(difference_type _Off)
		{	
		*(_Mybase *)this += _Off;
		return (*this);
		}

	_Myiter operator+(difference_type _Off) const
		{	
		_Myiter _Tmp = *this;
		return (_Tmp += _Off);
		}

	_Myiter& operator-=(difference_type _Off)
		{	
		return (*this += -_Off);
		}

	_Myiter operator-(difference_type _Off) const
		{	
		_Myiter _Tmp = *this;
		return (_Tmp -= _Off);
		}

	difference_type operator-(const _Mybase& _Right) const
		{	
		return (*(_Mybase *)this - _Right);
		}

	reference operator[](difference_type _Off) const
		{	
		return (*(*this + _Off));
		}
	};

template<class _Ty,
	size_t _Size> inline
	typename _Array_iterator<_Ty, _Size>::_Unchecked_type
		_Unchecked(_Array_iterator<_Ty, _Size> _Iter)
	{	
	return (_Iter._Unchecked());
	}

template<class _Ty,
	size_t _Size> inline
	_Array_iterator<_Ty, _Size>&
		_Rechecked(_Array_iterator<_Ty, _Size>& _Iter,
			typename _Array_iterator<_Ty, _Size>
				::_Unchecked_type _Right)
	{	
	return (_Iter._Rechecked(_Right));
	}

template<class _Ty,
	size_t _Size> inline
	_Array_iterator<_Ty, _Size> operator+(
		typename _Array_iterator<_Ty, _Size>::difference_type _Off,
		_Array_iterator<_Ty, _Size> _Next)
	{	
	return (_Next += _Off);
	}

		
template<class _RanIt>
	class move_iterator
	{	
public:
	typedef move_iterator<_RanIt> _Myt;
	typedef typename iterator_traits<_RanIt>::iterator_category
		iterator_category;
	typedef typename iterator_traits<_RanIt>::value_type
		value_type;
	typedef typename iterator_traits<_RanIt>::difference_type
		difference_type;
	typedef _RanIt pointer;
	typedef typename iterator_traits<_RanIt>::reference _Ref0;
	typedef conditional_t<is_reference<_Ref0>::value,
		remove_reference_t<_Ref0>&&, _Ref0> reference;
	typedef _RanIt iterator_type;

	move_iterator()
		: current()
		{	
		}

	explicit move_iterator(iterator_type _Right)
		: current(_Right)
		{	
		}

	template<class _RanIt2>
		move_iterator(const move_iterator<_RanIt2>& _Right)
		: current(_Right.base())
		{	
		}

	template<class _RanIt2>
		_Myt& operator=(const move_iterator<_RanIt2>& _Right)
		{	
		current = _Right.base();
		return (*this);
		}

	_RanIt base() const
		{	
		return (current);
		}

	reference operator*() const
		{	
		return (static_cast<reference>(*current));
		}

	pointer operator->() const
		{	
		return (current);
		}

	_Myt& operator++()
		{	
		++current;
		return (*this);
		}

	_Myt operator++(int)
		{	
		_Myt _Tmp = *this;
		++current;
		return (_Tmp);
		}

	_Myt& operator--()
		{	
		--current;
		return (*this);
		}

	_Myt operator--(int)
		{	
		_Myt _Tmp = *this;
		--current;
		return (_Tmp);
		}

	template<class _RanIt2>
		bool _Equal(const move_iterator<_RanIt2>& _Right) const
		{	
		return (current == _Right.base());
		}



	_Myt& operator+=(difference_type _Off)
		{	
		current += _Off;
		return (*this);
		}

	_Myt operator+(difference_type _Off) const
		{	
		return (_Myt(current + _Off));
		}

	_Myt& operator-=(difference_type _Off)
		{	
		current -= _Off;
		return (*this);
		}

	_Myt operator-(difference_type _Off) const
		{	
		return (_Myt(current - _Off));
		}

	reference operator[](difference_type _Off) const
		{	
		return (::std:: move(current[_Off]));
		}

	template<class _RanIt2>
		bool _Less(const move_iterator<_RanIt2>& _Right) const
		{	
		return (current < _Right.base());
		}

	difference_type operator-(const _Myt& _Right) const
		{	
		return (current - _Right.base());
		}

protected:
	iterator_type current;	
	};

template<class _RanIt>
	struct _Is_checked_helper<move_iterator<_RanIt> >
		: public _Is_checked_helper<_RanIt>
	{	
	};

		
template<class _RanIt,
	class _Diff> inline
	move_iterator<_RanIt>
		operator+(_Diff _Off,
		const move_iterator<_RanIt>& _Right)
	{	
	return (_Right + _Off);
	}

template<class _RanIt1,
	class _RanIt2>
	auto inline operator-(
		move_iterator<_RanIt1>& _Left,
		const move_iterator<_RanIt2>& _Right)
			-> decltype(_Left.base() - _Right.base())
	{	
	return (_Left.base() - _Right.base());
	}

template<class _RanIt1,
	class _RanIt2> inline
	bool operator==(
		const move_iterator<_RanIt1>& _Left,
		const move_iterator<_RanIt2>& _Right)
	{	
	return (_Left._Equal(_Right));
	}

template<class _RanIt1,
	class _RanIt2> inline
	bool operator!=(
		const move_iterator<_RanIt1>& _Left,
		const move_iterator<_RanIt2>& _Right)
	{	
	return (!(_Left == _Right));
	}

template<class _RanIt1,
	class _RanIt2> inline
	bool operator<(
		const move_iterator<_RanIt1>& _Left,
		const move_iterator<_RanIt2>& _Right)
	{	
	return (_Left._Less(_Right));
	}

template<class _RanIt1,
	class _RanIt2> inline
	bool operator>(
		const move_iterator<_RanIt1>& _Left,
		const move_iterator<_RanIt2>& _Right)
	{	
	return (_Right < _Left);
	}

template<class _RanIt1,
	class _RanIt2> inline
	bool operator<=(
		const move_iterator<_RanIt1>& _Left,
		const move_iterator<_RanIt2>& _Right)
	{	
	return (!(_Right < _Left));
	}

template<class _RanIt1,
	class _RanIt2> inline
	bool operator>=(
		const move_iterator<_RanIt1>& _Left,
		const move_iterator<_RanIt2>& _Right)
	{	
	return (!(_Left < _Right));
	}

		
template<class _RanIt> inline
	move_iterator<_RanIt> make_move_iterator(_RanIt _Iter)
	{	
	return (move_iterator<_RanIt>(_Iter));
	}

		
template<class _Traits>
	struct _Char_traits_eq
	{
	typedef typename _Traits::char_type _Elem;

	bool operator()(_Elem _Left, _Elem _Right) const
		{
		return (_Traits::eq(_Left, _Right));
		}
	};

		
template<class _Traits>
	struct _Char_traits_lt
	{
	typedef typename _Traits::char_type _Elem;

	bool operator()(_Elem _Left, _Elem _Right) const
		{
		return (_Traits::lt(_Left, _Right));
		}
	};

		
template<class _InIt,
	class _OutIt> inline
	_OutIt _Copy_memmove(_InIt _First, _InIt _Last,
		_OutIt _Dest)
	{	
	const char * const _First_ch = reinterpret_cast<const char *>(_First);
	const char * const _Last_ch = reinterpret_cast<const char *>(_Last);
	char * const _Dest_ch = reinterpret_cast<char *>(_Dest);
	const size_t _Count = _Last_ch - _First_ch;
	:: memmove(_Dest_ch, _First_ch, _Count);
	return (reinterpret_cast<_OutIt>(_Dest_ch + _Count));
	}

template<class _InIt,
	class _OutIt> inline
	_OutIt _Copy_unchecked1(_InIt _First, _InIt _Last,
		_OutIt _Dest, _General_ptr_iterator_tag)
	{	
	for (; _First != _Last; ++_Dest, (void)++_First)
		*_Dest = *_First;
	return (_Dest);
	}

template<class _InIt,
	class _OutIt> inline
	_OutIt _Copy_unchecked1(_InIt _First, _InIt _Last,
		_OutIt _Dest, _Trivially_copyable_ptr_iterator_tag)
	{	
	return (_Copy_memmove(_First, _Last, _Dest));
	}

template<class _InIt,
	class _OutIt> inline
	_OutIt _Copy_unchecked(_InIt _First, _InIt _Last,
		_OutIt _Dest)
	{	
		
	return (_Copy_unchecked1(_First, _Last,
		_Dest, _Ptr_copy_cat(_First, _Dest)));
	}

template<class _InIt,
	class _OutIt> inline
	_OutIt _Copy_no_deprecate1(_InIt _First, _InIt _Last,
		_OutIt _Dest, input_iterator_tag, _Any_tag)
	{	
	return (_Rechecked(_Dest,
		_Copy_unchecked(_First, _Last, _Unchecked_idl0(_Dest))));
	}

template<class _InIt,
	class _OutIt> inline
	_OutIt _Copy_no_deprecate1(_InIt _First, _InIt _Last,
		_OutIt _Dest, random_access_iterator_tag, random_access_iterator_tag)
	{	
	;
	return (_Rechecked(_Dest,
		_Copy_unchecked(_First, _Last, _Unchecked(_Dest))));
	}

template<class _InIt,
	class _OutIt> inline
	_OutIt _Copy_no_deprecate(_InIt _First, _InIt _Last,
		_OutIt _Dest)
	{	
	;
	return (_Copy_no_deprecate1(_Unchecked(_First), _Unchecked(_Last),
		_Dest, _Iter_cat_t<_InIt>(), _Iter_cat_t<_OutIt>()));
	}

template<class _InIt,
	class _OutIt> inline
	_OutIt copy(_InIt _First, _InIt _Last,
		_OutIt _Dest)
	{	
	struct _Unchecked_iterators { static void  _Deprecate(false_type) { } static void _Deprecate(true_type) { } }; (_Unchecked_iterators::_Deprecate(_Is_checked(_Dest)));
	return (_Copy_no_deprecate(_First, _Last, _Dest));
	}

 












		
template<class _InIt,
	class _Diff,
	class _OutIt> inline
	_OutIt _Copy_n_unchecked2(_InIt _First, _Diff _Count,
		_OutIt _Dest, input_iterator_tag)
	{	
	if (0 < _Count)
		{
		*_Dest = *_First;
		while (0 < --_Count)
			*++_Dest = *++_First;
		return (++_Dest);
		}

	return (_Dest);
	}

template<class _InIt,
	class _Diff,
	class _OutIt> inline
	_OutIt _Copy_n_unchecked2(_InIt _First, _Diff _Count,
		_OutIt _Dest, forward_iterator_tag)
	{	
	for (; 0 < _Count; --_Count, (void)++_Dest, ++_First)
		*_Dest = *_First;
	return (_Dest);
	}

template<class _InIt,
	class _Diff,
	class _OutIt> inline
	_OutIt _Copy_n_unchecked1(_InIt _First, _Diff _Count,
		_OutIt _Dest, _General_ptr_iterator_tag)
	{	
		
		
	return (_Copy_n_unchecked2(_First, _Count,
		_Dest, _Iter_cat_t<_InIt>()));
	}

template<class _InIt,
	class _Diff,
	class _OutIt> inline
	_OutIt _Copy_n_unchecked1(_InIt _First, _Diff _Count,
		_OutIt _Dest, _Trivially_copyable_ptr_iterator_tag)
	{	
	if (0 < _Count)
		return (_Copy_memmove(_First, _First + _Count, _Dest));
	return (_Dest);
	}

template<class _InIt,
	class _Diff,
	class _OutIt> inline
	_OutIt _Copy_n_unchecked(_InIt _First, _Diff _Count,
		_OutIt _Dest)
	{	
	return (_Copy_n_unchecked1(_First, _Count,
		_Dest, _Ptr_copy_cat(_First, _Dest)));
	}

template<class _InIt,
	class _Diff,
	class _OutIt> inline
	_OutIt copy_n(_InIt _First, _Diff _Count,
		_OutIt _Dest)
	{	
		
	struct _Unchecked_iterators { static void  _Deprecate(false_type) { } static void _Deprecate(true_type) { } }; (_Unchecked_iterators::_Deprecate(_Is_checked(_Dest)));
	return (_Rechecked(_Dest,
		_Copy_n_unchecked(_Unchecked_n(_First, _Count), _Count, _Unchecked_n(_Dest, _Count))));
	}

 







































		
template<class _BidIt1,
	class _BidIt2> inline
	_BidIt2 _Copy_backward_memmove(_BidIt1 _First, _BidIt1 _Last,
		_BidIt2 _Dest)
	{	
	const char * const _First_ch = reinterpret_cast<const char *>(_First);
	const char * const _Last_ch = reinterpret_cast<const char *>(_Last);
	char * const _Dest_ch = reinterpret_cast<char *>(_Dest);
	const size_t _Count = _Last_ch - _First_ch;
	return (static_cast<_BidIt2>(
		:: memmove(_Dest_ch - _Count, _First_ch, _Count)));
	}

template<class _BidIt1,
	class _BidIt2> inline
	_BidIt2 _Copy_backward_unchecked1(_BidIt1 _First, _BidIt1 _Last,
		_BidIt2 _Dest, _General_ptr_iterator_tag)
	{	
	while (_First != _Last)
		*--_Dest = *--_Last;
	return (_Dest);
	}

template<class _BidIt1,
	class _BidIt2> inline
	_BidIt2 _Copy_backward_unchecked1(_BidIt1 _First, _BidIt1 _Last,
		_BidIt2 _Dest, _Trivially_copyable_ptr_iterator_tag)
	{	
	return (_Copy_backward_memmove(_First, _Last, _Dest));
	}

template<class _BidIt1,
	class _BidIt2> inline
	_BidIt2 _Copy_backward_unchecked(_BidIt1 _First, _BidIt1 _Last,
		_BidIt2 _Dest)
	{	
	return (_Copy_backward_unchecked1(_First, _Last,
		_Dest, _Ptr_copy_cat(_First, _Dest)));
	}

template<class _BidIt1,
	class _BidIt2> inline
	_BidIt2 _Copy_backward1(_BidIt1 _First, _BidIt1 _Last,
		_BidIt2 _Dest, input_iterator_tag, _Any_tag)
	{	
	return (_Rechecked(_Dest,
		_Copy_backward_unchecked(_First, _Last, _Unchecked_idl0(_Dest))));
	}

template<class _BidIt1,
	class _BidIt2> inline
	_BidIt2 _Copy_backward1(_BidIt1 _First, _BidIt1 _Last,
		_BidIt2 _Dest, random_access_iterator_tag, random_access_iterator_tag)
	{	
	;
	return (_Rechecked(_Dest,
		_Copy_backward_unchecked(_First, _Last, _Unchecked(_Dest))));
	}

template<class _BidIt1,
	class _BidIt2> inline
	_BidIt2 copy_backward(_BidIt1 _First, _BidIt1 _Last,
		_BidIt2 _Dest)
	{	
	struct _Unchecked_iterators { static void  _Deprecate(false_type) { } static void _Deprecate(true_type) { } }; (_Unchecked_iterators::_Deprecate(_Is_checked(_Dest)));
	;
	return (_Copy_backward1(_Unchecked(_First), _Unchecked(_Last),
		_Dest, _Iter_cat_t<_BidIt1>(), _Iter_cat_t<_BidIt2>()));
	}

		
template<class _InIt,
	class _OutIt> inline
	_OutIt _Move_unchecked1(_InIt _First, _InIt _Last,
		_OutIt _Dest, _General_ptr_iterator_tag)
	{	
	for (; _First != _Last; ++_Dest, (void)++_First)
		*_Dest = ::std:: move(*_First);
	return (_Dest);
	}

template<class _InIt,
	class _OutIt> inline
	_OutIt _Move_unchecked1(_InIt _First, _InIt _Last,
		_OutIt _Dest, _Trivially_copyable_ptr_iterator_tag)
	{	
	return (_Copy_memmove(_First, _Last, _Dest));
	}

template<class _InIt,
	class _OutIt> inline
	_OutIt _Move_unchecked(_InIt _First, _InIt _Last,
		_OutIt _Dest)
	{	
	return (_Move_unchecked1(_First, _Last,
		_Dest, _Ptr_move_cat(_First, _Dest)));
	}

template<class _InIt,
	class _OutIt> inline
	_OutIt _Move_no_deprecate1(_InIt _First, _InIt _Last,
		_OutIt _Dest, input_iterator_tag, _Any_tag)
	{	
	return (_Rechecked(_Dest,
		_Move_unchecked(_First, _Last, _Unchecked_idl0(_Dest))));
	}

template<class _InIt,
	class _OutIt> inline
	_OutIt _Move_no_deprecate1(_InIt _First, _InIt _Last,
		_OutIt _Dest, random_access_iterator_tag, random_access_iterator_tag)
	{	
	;
	return (_Rechecked(_Dest,
		_Move_unchecked(_First, _Last, _Unchecked(_Dest))));
	}

template<class _InIt,
	class _OutIt> inline
	_OutIt _Move_no_deprecate(_InIt _First, _InIt _Last,
		_OutIt _Dest)
	{	
	;
	return (_Move_no_deprecate1(_Unchecked(_First), _Unchecked(_Last),
		_Dest, _Iter_cat_t<_InIt>(), _Iter_cat_t<_OutIt>()));
	}

template<class _InIt,
	class _OutIt> inline
	_OutIt move(_InIt _First, _InIt _Last,
		_OutIt _Dest)
	{	
	struct _Unchecked_iterators { static void  _Deprecate(false_type) { } static void _Deprecate(true_type) { } }; (_Unchecked_iterators::_Deprecate(_Is_checked(_Dest)));
	return (_Move_no_deprecate(_First, _Last, _Dest));
	}

 












		
template<class _BidIt1,
	class _BidIt2> inline
	_BidIt2 _Move_backward_unchecked1(_BidIt1 _First, _BidIt1 _Last,
		_BidIt2 _Dest, _General_ptr_iterator_tag)
	{	
	while (_First != _Last)
		*--_Dest = ::std:: move(*--_Last);
	return (_Dest);
	}

template<class _BidIt1,
	class _BidIt2> inline
	_BidIt2 _Move_backward_unchecked1(_BidIt1 _First, _BidIt1 _Last,
		_BidIt2 _Dest, _Trivially_copyable_ptr_iterator_tag)
	{	
	return (_Copy_backward_memmove(_First, _Last, _Dest));
	}

template<class _BidIt1,
	class _BidIt2> inline
	_BidIt2 _Move_backward_unchecked(_BidIt1 _First, _BidIt1 _Last,
		_BidIt2 _Dest)
	{	
	return (_Move_backward_unchecked1(_First, _Last,
		_Dest, _Ptr_move_cat(_First, _Dest)));
	}

template<class _BidIt1,
	class _BidIt2> inline
	_BidIt2 _Move_backward1(_BidIt1 _First, _BidIt1 _Last,
		_BidIt2 _Dest, input_iterator_tag, _Any_tag)
	{	
	return (_Rechecked(_Dest,
		_Move_backward_unchecked(_First, _Last, _Unchecked_idl0(_Dest))));
	}

template<class _BidIt1,
	class _BidIt2> inline
	_BidIt2 _Move_backward1(_BidIt1 _First, _BidIt1 _Last,
		_BidIt2 _Dest, random_access_iterator_tag, random_access_iterator_tag)
	{	
	;
	return (_Rechecked(_Dest,
		_Move_backward_unchecked(_First, _Last, _Unchecked(_Dest))));
	}

template<class _BidIt1,
	class _BidIt2> inline
	_BidIt2 move_backward(_BidIt1 _First, _BidIt1 _Last,
		_BidIt2 _Dest)
	{	
	struct _Unchecked_iterators { static void  _Deprecate(false_type) { } static void _Deprecate(true_type) { } }; (_Unchecked_iterators::_Deprecate(_Is_checked(_Dest)));
	;
	return (_Move_backward1(_Unchecked(_First), _Unchecked(_Last),
		_Dest, _Iter_cat_t<_BidIt1>(), _Iter_cat_t<_BidIt2>()));
	}

		
template<class _Ty>
	struct _Is_character
		: false_type
		{	
		};

template<>
	struct _Is_character<char>
		: true_type
		{	
		};

template<>
	struct _Is_character<signed char>
		: true_type
		{	
		};

template<>
	struct _Is_character<unsigned char>
		: true_type
		{	
		};

template<class _FwdIt,
	class _Ty>
	struct _Fill_memset_is_safe_helper
	{	
	typedef _Iter_value_t<_FwdIt> _Value_type;
	typedef typename conjunction<
		is_pointer<_FwdIt>,
		disjunction<
			conjunction<
				_Is_character<_Ty>,
				_Is_character<_Value_type>>,
			conjunction<
				is_same<bool, _Ty>,
				is_same<bool, _Value_type>>
		>>::type type;
	};

template<class _FwdIt,
	class _Ty> inline
	typename _Fill_memset_is_safe_helper<_FwdIt, _Ty>::type
	_Fill_memset_is_safe(const _FwdIt&, const _Ty&)
	{	
	return {};
	}

template<class _FwdIt,
	class _Ty> inline
	void _Fill_unchecked1(_FwdIt _First, _FwdIt _Last, const _Ty& _Val, false_type)
	{	
	for (; _First != _Last; ++_First)
		*_First = _Val;
	}

template<class _FwdIt,
	class _Ty> inline
	void _Fill_unchecked1(_FwdIt _First, _FwdIt _Last, const _Ty& _Val, true_type)
	{	
	:: memset(_First, _Val, _Last - _First);
	}

template<class _FwdIt,
	class _Ty> inline
	void _Fill_unchecked(_FwdIt _First, _FwdIt _Last, const _Ty& _Val)
	{	
	_Fill_unchecked1(_First, _Last, _Val, _Fill_memset_is_safe(_First, _Val));
	}

template<class _FwdIt,
	class _Ty> inline
	void fill(_FwdIt _First, _FwdIt _Last, const _Ty& _Val)
	{	
	;
	_Fill_unchecked(_Unchecked(_First), _Unchecked(_Last), _Val);
	}

		
template<class _OutIt,
	class _Diff,
	class _Ty> inline
	_OutIt _Fill_n_unchecked1(_OutIt _Dest, _Diff _Count, const _Ty& _Val, false_type)
	{	
	for (; 0 < _Count; --_Count, (void)++_Dest)
		*_Dest = _Val;
	return (_Dest);
	}

template<class _OutIt,
	class _Diff,
	class _Ty> inline
	_OutIt _Fill_n_unchecked1(_OutIt _Dest, _Diff _Count, const _Ty& _Val, true_type)
	{	
	if (0 < _Count)
		{
		:: memset(_Dest, _Val, _Count);
		return (_Dest + _Count);
		}

	return (_Dest);
	}

template<class _OutIt,
	class _Diff,
	class _Ty> inline
	_OutIt _Fill_n_unchecked(_OutIt _Dest, _Diff _Count, const _Ty& _Val)
	{	
		
	return (_Fill_n_unchecked1(_Dest, _Count, _Val, _Fill_memset_is_safe(_Dest, _Val)));
	}

template<class _OutIt,
	class _Diff,
	class _Ty> inline
	_OutIt fill_n(_OutIt _Dest, _Diff _Count, const _Ty& _Val)
	{	
	return (_Rechecked(_Dest,
		_Fill_n_unchecked(_Unchecked_n(_Dest, _Count), _Count, _Val)));
	}

		
template<class _Elem1,
	class _Elem2>
	struct _Value_equality_is_bitwise_equality
		: bool_constant<static_cast<_Elem1>(-1) == static_cast<_Elem2>(-1)>
	{	
		
		
		
	};

template<class _Elem1,
	class _Elem2,
	class _Pr>
	struct _Equal_memcmp_is_safe_helper
		: false_type
	{	
		
	};

template<class _Elem1,
	class _Elem2>
	struct _Equal_memcmp_is_safe_helper<_Elem1, _Elem2, equal_to<>>
		: conjunction<
			_Is_same_size<_Elem1, _Elem2>,
			is_integral<_Elem1>,
			is_integral<_Elem2>,
			negation<is_same<bool, _Elem1>>,
			negation<is_same<bool, _Elem2>>,
			negation<is_volatile<_Elem1>>,
			negation<is_volatile<_Elem2>>,
			
			
			_Value_equality_is_bitwise_equality<_Elem1, _Elem2>
		>::type
	{	
	};

template<class _Elem1,
	class _Elem2>
	struct _Equal_memcmp_is_safe_helper<_Elem1 *, _Elem2 *, equal_to<>>
		: is_same<remove_cv_t<_Elem1>, remove_cv_t<_Elem2>>::type
	{	
	};

template<class _Elem>
	struct _Equal_memcmp_is_safe_helper<_Elem, _Elem, _Char_traits_eq<char_traits<_Elem>>>
		: _Equal_memcmp_is_safe_helper<_Elem, _Elem, equal_to<>>::type
	{	
	};

template<class _Elem>
	struct _Equal_memcmp_is_safe_helper<_Elem, _Elem, equal_to<_Elem>>
		: _Equal_memcmp_is_safe_helper<_Elem, _Elem, equal_to<>>::type
	{	
		
	};

template<class _Iter1,
	class _Iter2,
	class _Pr> inline
	false_type _Equal_memcmp_is_safe(const _Iter1&, const _Iter2&, const _Pr&)
	{	
	return {};
	}

template<class _Obj1,
	class _Obj2,
	class _Pr> inline
	typename _Equal_memcmp_is_safe_helper<
		remove_const_t<_Obj1>,
		remove_const_t<_Obj2>,
		_Pr>::type
		_Equal_memcmp_is_safe(_Obj1 * const&, _Obj2 * const&, const _Pr&)
	{	
	return {};
	}

template<class _InIt1,
	class _InIt2,
	class _Pr> inline
	bool _Equal_unchecked1(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _Pr& _Pred, false_type)
	{	
	for (; _First1 != _Last1; ++_First1, (void)++_First2)
		if (!_Pred(*_First1, *_First2))
			return (false);
	return (true);
	}

template<class _InIt1,
	class _InIt2,
	class _Pr> inline
	bool _Equal_unchecked1(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _Pr&, true_type)
	{	
	const char * const _First1_ch = reinterpret_cast<const char *>(_First1);
	const char * const _First2_ch = reinterpret_cast<const char *>(_First2);
	const size_t _Count = reinterpret_cast<const char *>(_Last1) - _First1_ch;
	return (:: memcmp(_First1_ch, _First2_ch, _Count) == 0);
	}

template<class _InIt1,
	class _InIt2,
	class _Pr> inline
	bool _Equal_unchecked(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _Pr& _Pred)
	{	
	return (_Equal_unchecked1(_First1, _Last1, _First2, _Pred,
		_Equal_memcmp_is_safe(_First1, _First2, _Pred)));
	}

template<class _InIt1,
	class _InIt2,
	class _Pr> inline
	bool _Equal_no_deprecate1(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _Pr& _Pred, input_iterator_tag, input_iterator_tag)
	{	
	return (_Equal_unchecked(_First1, _Last1, _Unchecked_idl0(_First2), _Pred));
	}

template<class _InIt1,
	class _InIt2,
	class _Pr> inline
	bool _Equal_no_deprecate1(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _Pr& _Pred, random_access_iterator_tag, random_access_iterator_tag)
	{	
	;
	return (_Equal_unchecked(_First1, _Last1, _Unchecked(_First2), _Pred));
	}

template<class _InIt1,
	class _InIt2,
	class _Pr> inline
	bool _Equal_no_deprecate(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _Pr& _Pred)
	{	
	;
	;
	return (_Equal_no_deprecate1(_Unchecked(_First1), _Unchecked(_Last1),
		_First2, _Pred, _Iter_cat_t<_InIt1>(), _Iter_cat_t<_InIt2>()));
	}

template<class _InIt1,
	class _InIt2,
	class _Pr> inline
	bool equal(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _Pr _Pred)
	{	
	struct _Unchecked_iterators { static void  _Deprecate(false_type) { } static void _Deprecate(true_type) { } }; (_Unchecked_iterators::_Deprecate(_Is_checked(_First2)));
	return (_Equal_no_deprecate(_First1, _Last1, _First2, _Pred));
	}

 













		
template<class _InIt1,
	class _InIt2> inline
	bool equal(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2)
	{	
	return (::std:: equal(_First1, _Last1, _First2,
		equal_to<>()));
	}

 











		
template<class _InIt1,
	class _InIt2,
	class _Pr> inline
	bool _Equal_unchecked(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _InIt2 _Last2, _Pr& _Pred,
			input_iterator_tag, input_iterator_tag)
	{	
		
	;
	for (; _First1 != _Last1 && _First2 != _Last2; ++_First1, (void)++_First2)
		if (!_Pred(*_First1, *_First2))
			return (false);
	return (_First1 == _Last1 && _First2 == _Last2);
	}

template<class _InIt1,
	class _InIt2,
	class _Pr> inline
	bool _Equal_unchecked(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _InIt2 _Last2, _Pr& _Pred,
			random_access_iterator_tag, random_access_iterator_tag)
	{	
		
	if (_Last1 - _First1 != _Last2 - _First2)
		return (false);
	;
	return (_Equal_unchecked(_First1, _Last1, _First2, _Pred));
	}

template<class _InIt1,
	class _InIt2,
	class _Pr> inline
	bool equal(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _InIt2 _Last2, _Pr _Pred)
	{	
	;
	;
	return (_Equal_unchecked(_Unchecked(_First1), _Unchecked(_Last1),
		_Unchecked(_First2), _Unchecked(_Last2), _Pred,
			_Iter_cat_t<_InIt1>(), _Iter_cat_t<_InIt2>()));
	}

		
template<class _InIt1,
	class _InIt2> inline
	bool equal(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _InIt2 _Last2)
	{	
	return (::std:: equal(_First1, _Last1, _First2, _Last2,
		equal_to<>()));
	}

		
template<class _Elem1,
	class _Elem2,
	class _FTy>
	struct _Lex_compare_check_element_types_helper
		: conjunction<
			_Is_character<_Elem1>,
			_Is_character<_Elem2>,
			_Is_character<_FTy>,
			is_unsigned<_FTy>
		>::type
	{	
	};

template<class _Elem1,
	class _Elem2>
	struct _Lex_compare_check_element_types_helper<_Elem1, _Elem2, void>
		: conjunction<
			_Is_character<_Elem1>,
			_Is_character<_Elem2>,
			is_unsigned<_Elem1>,
			is_unsigned<_Elem2>
		>::type
	{	
	};

template<class _Memcmp_pr>
	struct _Lex_compare_optimize
	{	
	};

template<class _Memcmp_pr,
	class _Obj1,
	class _Obj2,
	class _FTy>
	using _Lex_compare_check_element_types = _Lex_compare_optimize<conditional_t<
		_Lex_compare_check_element_types_helper<remove_const_t<_Obj1>, remove_const_t<_Obj2>, _FTy>::value,
		_Memcmp_pr, void>>;	

template<class _InIt1,
	class _InIt2,
	class _Pr> inline
	_Lex_compare_optimize<void> _Lex_compare_memcmp_classify(const _InIt1&, const _InIt2&, const _Pr&)
	{	
		
	return {};
	}

template<class _Obj1,
	class _Obj2,
	class _FTy> inline
	_Lex_compare_check_element_types<less<int>, _Obj1, _Obj2, _FTy>
		_Lex_compare_memcmp_classify(_Obj1 * const&, _Obj2 * const&, const less<_FTy>&)
	{	
	return {};
	}

template<class _InIt1,
	class _InIt2,
	class _Pr> inline
	bool _Lex_compare_unchecked1(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _InIt2 _Last2, _Pr& _Pred, _Lex_compare_optimize<void>)
	{	
	for (; _First1 != _Last1 && _First2 != _Last2; ++_First1, (void)++_First2)
		{	
		if (_Pred(*_First1, *_First2))
			return (true);
		else if (_Pred(*_First2, *_First1))
			return (false);
		}

	return (_First1 == _Last1 && _First2 != _Last2);
	}

template<class _InIt1,
	class _InIt2,
	class _Pr,
	class _Memcmp_pr> inline
	bool _Lex_compare_unchecked1(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _InIt2 _Last2, _Pr&, _Lex_compare_optimize<_Memcmp_pr>)
	{	
	const size_t _Num1 = _Last1 - _First1;
	const size_t _Num2 = _Last2 - _First2;
	const int _Ans = :: memcmp(_First1, _First2, _Num1 < _Num2 ? _Num1 : _Num2);
	return (_Memcmp_pr{}(_Ans, 0) || _Ans == 0 && _Num1 < _Num2);
	}

template<class _InIt1,
	class _InIt2,
	class _Pr> inline
	bool _Lex_compare_unchecked(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _InIt2 _Last2, _Pr& _Pred)
	{	
	return (_Lex_compare_unchecked1(_First1, _Last1, _First2, _Last2, _Pred,
		_Lex_compare_memcmp_classify(_First1, _First2, _Pred)));
	}

template<class _InIt1,
	class _InIt2,
	class _Pr> inline
	bool lexicographical_compare(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _InIt2 _Last2, _Pr _Pred)
	{	
	;
	;
	;
	return (_Lex_compare_unchecked(_Unchecked(_First1), _Unchecked(_Last1),
		_Unchecked(_First2), _Unchecked(_Last2), _Pred));
	}

		
template<class _InIt1,
	class _InIt2> inline
	bool lexicographical_compare(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _InIt2 _Last2)
	{	
	return (::std:: lexicographical_compare(_First1, _Last1,
		_First2, _Last2, less<>()));
	}

		
template<class _Ty> inline
	bool _Within_limits(const _Ty& _Val, true_type, true_type, _Any_tag)
	{	
	return ((-128) <= _Val && _Val <= 127);
	}

template<class _Ty> inline
	bool _Within_limits(const _Ty& _Val, true_type, false_type, true_type)
	{	
	return (_Val <= 127 || static_cast<_Ty>((-128)) <= _Val);
	}

template<class _Ty> inline
	bool _Within_limits(const _Ty& _Val, true_type, false_type, false_type)
	{	
	return (_Val <= 127);
	}

template<class _Ty> inline
	bool _Within_limits(const _Ty& _Val, false_type, true_type, _Any_tag)
	{	
	return (0 <= _Val && _Val <= 0xff);
	}

template<class _Ty> inline
	bool _Within_limits(const _Ty& _Val, false_type, false_type, _Any_tag)
	{	
	return (_Val <= 0xff);
	}

template<class _InIt,
	class _Ty> inline
	bool _Within_limits(_InIt, const _Ty& _Val)
	{	
	typedef typename remove_pointer<_InIt>::type _Elem;
	return (_Within_limits(_Val, is_signed<_Elem>(), is_signed<_Ty>(),
		integral_constant<bool, -1 == static_cast<_Ty>(-1)>()));
	}

template<class _InIt> inline
	bool _Within_limits(_InIt, const bool&)
	{	
	return (true);
	}

template<class _InIt,
	class _Ty> inline
	_InIt _Find_unchecked1(_InIt _First, _InIt _Last, const _Ty& _Val, true_type)
	{	
	if (!_Within_limits(_First, _Val))
		return (_Last);
	_First = static_cast<_InIt>(:: memchr(
		_First, static_cast<unsigned char>(_Val), _Last - _First));
	return (_First ? _First : _Last);
	}

template<class _InIt,
	class _Ty> inline
	_InIt _Find_unchecked1(_InIt _First, _InIt _Last, const _Ty& _Val, false_type)
	{	
	for (; _First != _Last; ++_First)
		if (*_First == _Val)
			break;
	return (_First);
	}

template<class _InIt,
	class _Ty> inline
	_InIt _Find_unchecked(_InIt _First, _InIt _Last, const _Ty& _Val)
	{	
	
	typedef integral_constant<bool,
		(is_same<_InIt, char *>::value
		|| is_same<_InIt, signed char *>::value
		|| is_same<_InIt, unsigned char *>::value
		|| is_same<_InIt, const char *>::value
		|| is_same<_InIt, const signed char *>::value
		|| is_same<_InIt, const unsigned char *>::value)
		&& is_integral<_Ty>::value
	> _Memchr_opt;
	return (_Find_unchecked1(_First, _Last, _Val, _Memchr_opt()));
	}

template<class _InIt,
	class _Ty> inline
	_InIt find(_InIt _First, _InIt _Last, const _Ty& _Val)
	{	
	;
	return (_Rechecked(_First,
		_Find_unchecked(_Unchecked(_First), _Unchecked(_Last), _Val)));
	}

		
template<class _InIt,
	class _Ty,
	class _Pr> inline
	_InIt _Find_pr(_InIt _First, _InIt _Last, const _Ty& _Val, _Pr& _Pred)
	{	
	for (; _First != _Last; ++_First)
		if (_Pred(*_First, _Val))
			break;
	return (_First);
	}

		
template<class _InIt,
	class _Ty> inline
	_Iter_diff_t<_InIt>
		_Count_unchecked(_InIt _First, _InIt _Last, const _Ty& _Val)
	{	
	_Iter_diff_t<_InIt> _Count = 0;

	for (; _First != _Last; ++_First)
		if (*_First == _Val)
			++_Count;
	return (_Count);
	}

template<class _InIt,
	class _Ty> inline
	_Iter_diff_t<_InIt>
		count(_InIt _First, _InIt _Last, const _Ty& _Val)
	{	
	;
	return (_Count_unchecked(_Unchecked(_First), _Unchecked(_Last), _Val));
	}

		
template<class _InIt,
	class _Ty,
	class _Pr> inline
	_Iter_diff_t<_InIt>
		_Count_pr(_InIt _First, _InIt _Last, const _Ty& _Val, _Pr& _Pred)
	{	
	_Iter_diff_t<_InIt> _Count = 0;

	for (; _First != _Last; ++_First)
		if (_Pred(*_First, _Val))
			++_Count;
	return (_Count);
	}

		
template<class _FwdIt1,
	class _FwdIt2,
	class _Pr> inline
	void _Trim_matching_suffixes(_FwdIt1&, _FwdIt2&, _Pr&,
		forward_iterator_tag, forward_iterator_tag)
	{	
	}

template<class _FwdIt1,
	class _FwdIt2,
	class _Pr> inline
	void _Trim_matching_suffixes(_FwdIt1& _Last1, _FwdIt2& _Last2, _Pr& _Pred,
		bidirectional_iterator_tag, bidirectional_iterator_tag)
	{	
	
	while (_Pred(*--_Last1, *--_Last2))
		;	
	++_Last1;
	++_Last2;
	}

		
template<class _FwdIt1,
	class _FwdIt2,
	class _Pr> inline
	bool _Check_match_counts(_FwdIt1 _First1, _FwdIt1 _Last1,
		_FwdIt2 _First2, _FwdIt2 _Last2, _Pr& _Pred)
	{	
	_Trim_matching_suffixes(_Last1, _Last2, _Pred,
		_Iter_cat_t<_FwdIt1>(), _Iter_cat_t<_FwdIt2>());
	for (_FwdIt1 _Next1 = _First1; _Next1 != _Last1; ++_Next1)
		if (_Next1 == _Find_pr(_First1, _Next1, *_Next1, _Pred))
			{	
			_Iter_diff_t<_FwdIt2> _Count2 = _Count_pr(_First2, _Last2, *_Next1, _Pred);
			if (_Count2 == 0)
				return (false);	
			_FwdIt1 _Skip1 = ::std:: next(_Next1);
			_Iter_diff_t<_FwdIt1> _Count1 = _Count_pr(_Skip1, _Last1, *_Next1, _Pred) + 1;
			if (_Count2 != _Count1)
				return (false);	
			}

	return (true);
	}

		
template<class _FwdIt1,
	class _FwdIt2,
	class _Pr> inline
	bool _Is_permutation_unchecked(_FwdIt1 _First1, _FwdIt1 _Last1,
		_FwdIt2 _First2, _Pr& _Pred)
	{	
	for (; _First1 != _Last1; ++_First1, (void)++_First2)
		if (!_Pred(*_First1, *_First2))
			{	
			_FwdIt2 _Last2 = ::std:: next(_First2,
				::std:: distance(_First1, _Last1));
			return (_Check_match_counts(_First1, _Last1,
				_First2, _Last2, _Pred));
			}

	return (true);
	}

template<class _FwdIt1,
	class _FwdIt2,
	class _Pr> inline
	bool _Is_permutation_no_deprecate1(_FwdIt1 _First1, _FwdIt1 _Last1,
		_FwdIt2 _First2, _Pr& _Pred, forward_iterator_tag, forward_iterator_tag)
	{	
	return (_Is_permutation_unchecked(_First1, _Last1, _Unchecked_idl0(_First2), _Pred));
	}

template<class _FwdIt1,
	class _FwdIt2,
	class _Pr> inline
	bool _Is_permutation_no_deprecate1(_FwdIt1 _First1, _FwdIt1 _Last1,
		_FwdIt2 _First2, _Pr& _Pred, random_access_iterator_tag, random_access_iterator_tag)
	{	
	;
	return (_Is_permutation_unchecked(_First1, _Last1, _Unchecked(_First2), _Pred));
	}

template<class _FwdIt1,
	class _FwdIt2,
	class _Pr> inline
	bool _Is_permutation_no_deprecate(_FwdIt1 _First1, _FwdIt1 _Last1,
		_FwdIt2 _First2, _Pr& _Pred)
	{	
	;
	;
	return (_Is_permutation_no_deprecate1(_Unchecked(_First1), _Unchecked(_Last1),
		_First2, _Pred, _Iter_cat_t<_FwdIt1>(), _Iter_cat_t<_FwdIt2>()));
	}

template<class _FwdIt1,
	class _FwdIt2,
	class _Pr> inline
	bool is_permutation(_FwdIt1 _First1, _FwdIt1 _Last1,
		_FwdIt2 _First2, _Pr _Pred)
	{	
	struct _Unchecked_iterators { static void  _Deprecate(false_type) { } static void _Deprecate(true_type) { } }; (_Unchecked_iterators::_Deprecate(_Is_checked(_First2)));
	return (_Is_permutation_no_deprecate(_First1, _Last1, _First2, _Pred));
	}

 













		
template<class _FwdIt1,
	class _FwdIt2> inline
	bool is_permutation(_FwdIt1 _First1, _FwdIt1 _Last1,
		_FwdIt2 _First2)
	{	
	return (::std:: is_permutation(_First1, _Last1,
		_First2, equal_to<>()));
	}


 










		
template<class _FwdIt1,
	class _FwdIt2,
	class _Pr> inline
	bool _Is_permutation_unchecked(_FwdIt1 _First1, _FwdIt1 _Last1,
		_FwdIt2 _First2, _FwdIt2 _Last2, _Pr& _Pred,
		forward_iterator_tag, forward_iterator_tag)
	{	
		
	;
	for (; _First1 != _Last1 && _First2 != _Last2; ++_First1, (void)++_First2)
		if (!_Pred(*_First1, *_First2))
			{	
			if (::std:: distance(_First1, _Last1)
				!= ::std:: distance(_First2, _Last2))
				return (false);	
			else
				return (_Check_match_counts(_First1, _Last1,
					_First2, _Last2, _Pred));
			}

	return (_First1 == _Last1 && _First2 == _Last2);
	}

template<class _FwdIt1,
	class _FwdIt2,
	class _Pr> inline
	bool _Is_permutation_unchecked(_FwdIt1 _First1, _FwdIt1 _Last1,
		_FwdIt2 _First2, _FwdIt2 _Last2, _Pr& _Pred,
		random_access_iterator_tag, random_access_iterator_tag)
	{	
		
	if (_Last1 - _First1 != _Last2 - _First2)
		return (false);
	;
	return (_Is_permutation_unchecked(_First1, _Last1, _First2, _Pred));
	}

template<class _FwdIt1,
	class _FwdIt2,
	class _Pr> inline
	bool is_permutation(_FwdIt1 _First1, _FwdIt1 _Last1,
		_FwdIt2 _First2, _FwdIt2 _Last2, _Pr _Pred)
	{	
		
	;
	;
	return (_Is_permutation_unchecked(_Unchecked(_First1), _Unchecked(_Last1),
		_Unchecked(_First2), _Unchecked(_Last2), _Pred,
		_Iter_cat_t<_FwdIt1>(), _Iter_cat_t<_FwdIt2>()));
	}

		
template<class _FwdIt1,
	class _FwdIt2> inline
	bool is_permutation(_FwdIt1 _First1, _FwdIt1 _Last1,
		_FwdIt2 _First2, _FwdIt2 _Last2)
	{	
	return (::std:: is_permutation(_First1, _Last1,
		_First2, _Last2, equal_to<>()));
	}

		
template<class _BidIt> inline
	void _Reverse_unchecked(_BidIt _First, _BidIt _Last)
	{	
	for (; _First != _Last && _First != --_Last; ++_First)
		::std:: iter_swap(_First, _Last);
	}

template<class _BidIt> inline
	void reverse(_BidIt _First, _BidIt _Last)
	{	
	;
	_Reverse_unchecked(_Unchecked(_First), _Unchecked(_Last));
	}

		
template<class _FwdIt> inline
	_FwdIt _Rotate_unchecked1(_FwdIt _First, _FwdIt _Mid, _FwdIt _Last,
		forward_iterator_tag)
	{	
	for (_FwdIt _Next = _Mid, _Res = _Last; ; )
		{	
		::std:: iter_swap(_First, _Next);
		if (++_First == _Mid)
			{	
			if (++_Next == _Last)
				return (_Res == _Last ? _Mid : _Res);
			else
				_Mid = _Next;	
			}
		else if (++_Next == _Last)
			{	
			if (_Res == _Last)
				_Res = _First;
			_Next = _Mid;
			}
		}
	}

template<class _BidIt> inline
	pair<_BidIt, _BidIt> _Reverse_until_sentinel_unchecked(
		_BidIt _First, _BidIt _Sentinel, _BidIt _Last)
	{	
	while (_First != _Sentinel && _Last != _Sentinel)
		::std:: iter_swap(_First++, --_Last);
	return (::std:: make_pair(_First, _Last));
	}

template<class _BidIt> inline
	_BidIt _Rotate_unchecked1(_BidIt _First, _BidIt _Mid, _BidIt _Last,
		bidirectional_iterator_tag)
	{	
	_Reverse_unchecked(_First, _Mid);
	_Reverse_unchecked(_Mid, _Last);
	pair<_BidIt, _BidIt> _Tmp = _Reverse_until_sentinel_unchecked(_First, _Mid, _Last);
	_Reverse_unchecked(_Tmp.first, _Tmp.second);
	return (_Mid != _Tmp.first ? _Tmp.first : _Tmp.second);
	}

template<class _RanIt> inline
	_RanIt _Rotate_unchecked1(_RanIt _First, _RanIt _Mid, _RanIt _Last,
		random_access_iterator_tag)
	{	
	_Reverse_unchecked(_First, _Mid);
	_Reverse_unchecked(_Mid, _Last);
	_Reverse_unchecked(_First, _Last);
	return (_First + (_Last - _Mid));
	}

template<class _FwdIt> inline
	_FwdIt _Rotate_unchecked(_FwdIt _First, _FwdIt _Mid, _FwdIt _Last)
	{	
	if (_First == _Mid)
		return (_Last);
	if (_Mid == _Last)
		return (_First);
	return (_Rotate_unchecked1(_First, _Mid, _Last, _Iter_cat_t<_FwdIt>()));
	}

template<class _FwdIt> inline
	_FwdIt rotate(_FwdIt _First, _FwdIt _Mid, _FwdIt _Last)
	{	
	;
	;
	return (_Rechecked(_First,
		_Rotate_unchecked(_Unchecked(_First), _Unchecked(_Mid),
		_Unchecked(_Last))));
	}

	
template<class _Diff,
	class _Urng>
	class _Rng_from_urng
	{	
public:
	typedef typename make_unsigned<_Diff>::type _Ty0;
	typedef typename _Urng::result_type _Ty1;

	typedef typename _If<sizeof (_Ty1) < sizeof (_Ty0),
		_Ty0, _Ty1>::type _Udiff;


	explicit _Rng_from_urng(_Urng& _Func)
		: _Ref(_Func), _Bits(8 * sizeof (_Udiff)), _Bmask(_Udiff(-1))
		{	
		for (; (_Urng::max)() - (_Urng::min)() < _Bmask; _Bmask >>= 1)
			--_Bits;
		}

	_Diff operator()(_Diff _Index)
		{	
		for (; ; )
			{	
			_Udiff _Ret = 0;	
			_Udiff _Mask = 0;	

			while (_Mask < _Udiff(_Index - 1))
				{	
				_Ret <<= _Bits - 1;	
				_Ret <<= 1;
				_Ret |= _Get_bits();
				_Mask <<= _Bits - 1;	
				_Mask <<= 1;
				_Mask |= _Bmask;
				}

			
			if (_Ret / _Index < _Mask / _Index
				|| _Mask % _Index == _Udiff(_Index - 1))
				return (_Ret % _Index);
			}
		}

	_Udiff _Get_all_bits()
		{	
		_Udiff _Ret = 0;

		for (size_t _Num = 0; _Num < 8 * sizeof (_Udiff);
			_Num += _Bits)
			{	
			_Ret <<= _Bits - 1;	
			_Ret <<= 1;
			_Ret |= _Get_bits();
			}

		return (_Ret);
		}

	_Rng_from_urng(const _Rng_from_urng&) = delete;
	_Rng_from_urng& operator=(const _Rng_from_urng&) = delete;

private:
	_Udiff _Get_bits()
		{	
		for (; ; )
			{	
			_Udiff _Val = _Ref() - (_Urng::min)();

			if (_Val <= _Bmask)
				return (_Val);
			}
		}

	_Urng& _Ref;	
	size_t _Bits;	
	_Udiff _Bmask;	
	};

		
template<class _Elem>
	class __declspec(dllimport) _Yarn
	{	
public:
	typedef _Yarn<_Elem> _Myt;

	 _Yarn()
		: _Myptr(0), _Nul(0)
		{	
		}

	 _Yarn(const _Myt& _Right)
		: _Myptr(0), _Nul(0)
		{	
		*this = _Right;
		}

	 _Yarn(const _Elem *_Right)
		: _Myptr(0), _Nul(0)
		{	
		*this = _Right;
		}

	_Myt&  operator=(const _Myt& _Right)
		{	
		return (*this = _Right._Myptr);
		}

	_Myt&  operator=(const _Elem *_Right)
		{	
		if (_Myptr != _Right)
			{	
			_Tidy();

			if (_Right != 0)
				{	
				const _Elem *_Ptr = _Right;
				while (*_Ptr != (_Elem)0)
					++_Ptr;
				size_t _Count = ((const char *)++_Ptr - (const char *)_Right);

 




				_Myptr = (_Elem *):: malloc(_Count);
 

				if (_Myptr != 0)
					:: memcpy(_Myptr, _Right, _Count);
				}
			}

		return (*this);
		}

	 ~_Yarn() noexcept
		{	
		_Tidy();
		}

	bool  empty() const
		{	
		return (_Myptr == 0);
		}

	const _Elem * c_str() const
		{	
		return (_Myptr != 0 ? _Myptr : &_Nul);
		}

	bool  _Empty() const
		{	
		return (_Myptr == 0);
		}

	const _Elem * _C_str() const
		{	
		return (_Myptr != 0 ? _Myptr : &_Nul);
		}

private:
	void  _Tidy()
		{	
		if (_Myptr != 0)

 



			:: free(_Myptr);
 

		_Myptr = 0;
		}

	_Elem *_Myptr;	
	_Elem _Nul;		
	};

	
template<class _Ty,
	class _Alloc>
	struct _Has_allocator_type
	{	
	template<class _Uty>
		static auto _Fn(int)
			-> is_convertible<_Alloc,
				typename _Uty::allocator_type>;
	template<class _Uty>
		static auto _Fn(_Wrap_int)
			-> false_type;

	typedef decltype(_Fn<_Ty>(0)) type;
	};

		
struct allocator_arg_t
	{	
	};

constexpr allocator_arg_t allocator_arg{};

[[noreturn]] __declspec(dllimport) void __cdecl _Xbad_alloc();
[[noreturn]] __declspec(dllimport) void __cdecl _Xinvalid_argument(  const char *);
[[noreturn]] __declspec(dllimport) void __cdecl _Xlength_error(  const char *);
[[noreturn]] __declspec(dllimport) void __cdecl _Xout_of_range(  const char *);
[[noreturn]] __declspec(dllimport) void __cdecl _Xoverflow_error(  const char *);
[[noreturn]] __declspec(dllimport) void __cdecl _Xruntime_error(  const char *);
}

namespace std {
		
template<class _Ty,
	class _Alloc>
	struct uses_allocator
		: _Has_allocator_type<_Ty, _Alloc>::type
	{	
	};

 
template<class _Ty,
	class _Alloc>
	constexpr bool uses_allocator_v = uses_allocator<_Ty, _Alloc>::value;
 
}	
 
 #pragma warning(pop)
 #pragma pack(pop)










 #pragma pack(push,8)
 #pragma warning(push,3)
 
 

 
  
 

 #pragma warning(disable: 4100)

namespace std {



 




 

 









		
inline
	__declspec(allocator) void *_Allocate(size_t _Count, size_t _Sz,
		bool _Try_aligned_allocation = true)
	{	
	void *_Ptr = 0;

	if (_Count == 0)
		return (_Ptr);

	
	if ((size_t)(-1) / _Sz < _Count)
		_Xbad_alloc();	
	const size_t _User_size = _Count * _Sz;

 
	if (_Try_aligned_allocation
		&& 4096 <= _User_size)
		{	
		static_assert(sizeof (void *) < 32,
			"Big allocations should at least match vector register size");
		const size_t _Block_size = (sizeof(void *) + 32 - 1) + _User_size;
		if (_Block_size <= _User_size)
			_Xbad_alloc();	
		const uintptr_t _Ptr_container =
			reinterpret_cast<uintptr_t>(::operator new(_Block_size));
		{ if (!(_Ptr_container != 0)) { ((void)0); ::_invalid_parameter_noinfo_noreturn(); } ; };
		_Ptr = reinterpret_cast<void *>((_Ptr_container + (sizeof(void *) + 32 - 1))
			& ~(32 - 1));
		static_cast<uintptr_t *>(_Ptr)[-1] = _Ptr_container;

 


		}
	else
 

		{	
		_Ptr = ::operator new(_User_size);
		{ if (!(_Ptr != 0)) { ((void)0); ::_invalid_parameter_noinfo_noreturn(); } ; };
		}
	return (_Ptr);
	}

		
inline
	void _Deallocate(void * _Ptr, size_t _Count, size_t _Sz)
	{	
 
	{ if (!(_Count <= (size_t)(-1) / _Sz)) { ((void)0); ::_invalid_parameter_noinfo_noreturn(); } ; };
	const size_t _User_size = _Count * _Sz;
	if (4096 <= _User_size)
		{	
		const uintptr_t _Ptr_user = reinterpret_cast<uintptr_t>(_Ptr);
		{ if (!((_Ptr_user & (32 - 1)) == 0)) { ((void)0); ::_invalid_parameter_noinfo_noreturn(); } ; };
		const uintptr_t _Ptr_ptr = _Ptr_user - sizeof(void *);
		const uintptr_t _Ptr_container =
			*reinterpret_cast<uintptr_t *>(_Ptr_ptr);

 







		
		{ if (!(_Ptr_container < _Ptr_user)) { ((void)0); ::_invalid_parameter_noinfo_noreturn(); } ; };

 




		{ if (!(sizeof(void *) <= _Ptr_user - _Ptr_container)) { ((void)0); ::_invalid_parameter_noinfo_noreturn(); } ; };
 

		{ if (!(_Ptr_user - _Ptr_container <= (sizeof(void *) + 32 - 1))) { ((void)0); ::_invalid_parameter_noinfo_noreturn(); } ; };

		_Ptr = reinterpret_cast<void *>(_Ptr_container);
		}
 

	::operator delete(_Ptr);
	}

		
template<class _Ty1,
	class _Ty2> inline
	void _Construct(_Ty1 *_Ptr, _Ty2&& _Val)
	{	
	void *_Vptr = _Ptr;
	::new (_Vptr) _Ty1(::std:: forward<_Ty2>(_Val));
	}

template<class _Ty1> inline
	void _Construct(_Ty1 *_Ptr)
	{	
	void *_Vptr = _Ptr;

	::new (_Vptr) _Ty1();
	}

		
template<class _Alty>
	struct _Is_simple_alloc
		: _Cat_base<is_same<typename _Alty::size_type, size_t>::value
		&& is_same<typename _Alty::difference_type, ptrdiff_t>::value
		&& is_same<typename _Alty::pointer,
			typename _Alty::value_type *>::value
		&& is_same<typename _Alty::const_pointer,
			const typename _Alty::value_type *>::value
		&& is_same<typename _Alty::reference,
			typename _Alty::value_type&>::value
		&& is_same<typename _Alty::const_reference,
			const typename _Alty::value_type&>::value>
	{	
	};

		
template<class _Value_type>
	struct _Simple_types
	{	
	typedef _Value_type value_type;
	typedef size_t size_type;
	typedef ptrdiff_t difference_type;
	typedef value_type *pointer;
	typedef const value_type *const_pointer;
	typedef value_type& reference;
	typedef const value_type& const_reference;
	};

		
template<class _Alty,
	class _Pointer>
	struct _Get_voidptr
	{	
	typedef typename _Alty::template rebind<void>::other _Alvoid;
	typedef typename _Alvoid::pointer type;
	};

template<class _Alty,
	class _Ty>
	struct _Get_voidptr<_Alty, _Ty *>
	{	
	typedef void *type;
	};

		
template<class _Ty>
	struct _Get_first_parameter;

template<template<class, class...> class _Ty,
	class _First,
	class... _Rest>
	struct _Get_first_parameter<_Ty<_First, _Rest...> >
	{	
	typedef _First type;
	};

		
template<class _Newfirst,
	class _Ty>
	struct _Replace_first_parameter;

template<class _Newfirst,
	template<class, class...> class _Ty,
	class _First,
	class... _Rest>
	struct _Replace_first_parameter<_Newfirst, _Ty<_First, _Rest...> >
	{	
	typedef _Ty<_Newfirst, _Rest...> type;
	};

		
template<class _Ty>
	struct _Get_element_type
	{ template<class _Uty> static auto _Fn(int) -> _Identity<typename _Uty::element_type>; template<class _Uty> static auto _Fn(_Wrap_int) -> _Identity<typename _Get_first_parameter<_Uty>::type>; typedef decltype(_Fn<_Ty>(0)) _Decltype; typedef typename _Decltype::type type; };

		
template<class _Ty>
	struct _Get_ptr_difference_type
	{ template<class _Uty> static auto _Fn(int) -> _Identity<typename _Uty::difference_type>; template<class _Uty> static auto _Fn(_Wrap_int) -> _Identity<ptrdiff_t>; typedef decltype(_Fn<_Ty>(0)) _Decltype; typedef typename _Decltype::type type; };

		
template<class _Ty,
	class _Other>
	struct _Get_rebind_type
	{ template<class _Uty> static auto _Fn(int) -> _Identity<typename _Uty::template rebind<_Other>::other>; template<class _Uty> static auto _Fn(_Wrap_int) -> _Identity<typename _Replace_first_parameter<_Other , _Uty>::type>; typedef decltype(_Fn<_Ty>(0)) _Decltype; typedef typename _Decltype::type type; };

		
template<class _Ty>
	struct pointer_traits
	{	
	typedef typename _Get_element_type<_Ty>::type element_type;
	typedef _Ty pointer;
	typedef typename _Get_ptr_difference_type<_Ty>::type difference_type;

	template<class _Other>
		using rebind = typename _Get_rebind_type<_Ty, _Other>::type;

	typedef typename _If<is_void<element_type>::value,
		char&,
		typename add_lvalue_reference<element_type>::type>::type _Reftype;

	static pointer pointer_to(_Reftype _Val)
		{	
		return (_Ty::pointer_to(_Val));
		}
	};

		
template<class _Ty>
	struct pointer_traits<_Ty *>
	{	
	typedef _Ty element_type;
	typedef _Ty *pointer;
	typedef ptrdiff_t difference_type;

	template<class _Other>
		using rebind = _Other *;

	typedef typename _If<is_void<_Ty>::value,
		char&,
		typename add_lvalue_reference<_Ty>::type>::type _Reftype;

	static pointer pointer_to(_Reftype _Val)
		{	
		return (::std:: addressof(_Val));
		}
	};


		
template<class _Ptrty> inline
	void _Destroy(_Ptrty _Ptr)
	{	
	typedef typename pointer_traits<_Ptrty>::element_type _Ty;
	_Ptr->~_Ty();
	}

		
template<class _Ptrty> inline
	auto _Const_cast(_Ptrty _Ptr)
	{	
	using _Elem = typename pointer_traits<_Ptrty>::element_type;
	using _Modifiable = remove_const_t<_Elem>;
	using _Dest = typename pointer_traits<_Ptrty>::template rebind<_Modifiable>;

	return (pointer_traits<_Dest>::pointer_to(const_cast<_Modifiable&>(*_Ptr)));
	}

template<class _Ty> inline
	auto _Const_cast(_Ty * _Ptr)
	{	
	return (const_cast<remove_const_t<_Ty> *>(_Ptr));
	}


		
template<class _Ty>
	struct _Get_pointer_type
	{ template<class _Uty> static auto _Fn(int) -> _Identity<typename _Uty::pointer>; template<class _Uty> static auto _Fn(_Wrap_int) -> _Identity<typename _Ty::value_type *>; typedef decltype(_Fn<_Ty>(0)) _Decltype; typedef typename _Decltype::type type; };

		
template<class _Ty>
	struct _Get_const_pointer_type
	{ template<class _Uty> static auto _Fn(int) -> _Identity<typename _Uty::const_pointer>; template<class _Uty> static auto _Fn(_Wrap_int) -> _Identity<typename pointer_traits<typename _Get_pointer_type<_Ty>::type> ::template rebind<const typename _Ty::value_type>>; typedef decltype(_Fn<_Ty>(0)) _Decltype; typedef typename _Decltype::type type; };

		
template<class _Ty>
	struct _Get_void_pointer_type
	{ template<class _Uty> static auto _Fn(int) -> _Identity<typename _Uty::void_pointer>; template<class _Uty> static auto _Fn(_Wrap_int) -> _Identity<typename pointer_traits<typename _Get_pointer_type<_Ty>::type> ::template rebind<void>>; typedef decltype(_Fn<_Ty>(0)) _Decltype; typedef typename _Decltype::type type; };

		
template<class _Ty>
	struct _Get_const_void_pointer_type
	{ template<class _Uty> static auto _Fn(int) -> _Identity<typename _Uty::const_void_pointer>; template<class _Uty> static auto _Fn(_Wrap_int) -> _Identity<typename pointer_traits<typename _Get_pointer_type<_Ty>::type> ::template rebind<const void>>; typedef decltype(_Fn<_Ty>(0)) _Decltype; typedef typename _Decltype::type type; };

		
template<class _Ty>
	struct _Get_difference_type
	{ template<class _Uty> static auto _Fn(int) -> _Identity<typename _Uty::difference_type>; template<class _Uty> static auto _Fn(_Wrap_int) -> _Identity<typename _Get_ptr_difference_type< typename _Get_pointer_type<_Ty>::type>::type>; typedef decltype(_Fn<_Ty>(0)) _Decltype; typedef typename _Decltype::type type; };

		
template<class _Ty>
	struct _Get_size_type
	{ template<class _Uty> static auto _Fn(int) -> _Identity<typename _Uty::size_type>; template<class _Uty> static auto _Fn(_Wrap_int) -> _Identity<typename make_unsigned< typename _Get_difference_type<_Ty>::type>::type>; typedef decltype(_Fn<_Ty>(0)) _Decltype; typedef typename _Decltype::type type; };

		
template<class _Ty>
	struct _Get_propagate_on_container_copy
	{ template<class _Uty> static auto _Fn(int) -> _Identity<typename _Uty::propagate_on_container_copy_assignment>; template<class _Uty> static auto _Fn(_Wrap_int) -> _Identity<false_type>; typedef decltype(_Fn<_Ty>(0)) _Decltype; typedef typename _Decltype::type type; };

		
template<class _Ty>
	struct _Get_propagate_on_container_move
	{ template<class _Uty> static auto _Fn(int) -> _Identity<typename _Uty::propagate_on_container_move_assignment>; template<class _Uty> static auto _Fn(_Wrap_int) -> _Identity<false_type>; typedef decltype(_Fn<_Ty>(0)) _Decltype; typedef typename _Decltype::type type; };

		
template<class _Ty>
	struct _Get_propagate_on_container_swap
	{ template<class _Uty> static auto _Fn(int) -> _Identity<typename _Uty::propagate_on_container_swap>; template<class _Uty> static auto _Fn(_Wrap_int) -> _Identity<false_type>; typedef decltype(_Fn<_Ty>(0)) _Decltype; typedef typename _Decltype::type type; };

		
template<class _Ty>
	struct _Get_is_always_equal
	{ template<class _Uty> static auto _Fn(int) -> _Identity<typename _Uty::is_always_equal>; template<class _Uty> static auto _Fn(_Wrap_int) -> _Identity<typename is_empty<_Ty>::type>; typedef decltype(_Fn<_Ty>(0)) _Decltype; typedef typename _Decltype::type type; };

		
template<class _Ty>
	class allocator;
template<class _Alloc>
	struct _Wrap_alloc;

template<class _Alloc>
	struct _Unwrap_alloc
	{	
	typedef _Alloc type;
	};

template<class _Alloc>
	struct _Unwrap_alloc<_Wrap_alloc<_Alloc>>
	{	
	typedef _Alloc type;
	};


		
template<class _Alloc>
	using _Unwrap_alloc_t = typename _Unwrap_alloc<_Alloc>::type;


		
template<class _Alloc,
	class = void>
	struct _Is_default_allocator
		: false_type
	{	
	};

template<class _Ty>
	struct _Is_default_allocator<allocator<_Ty>, typename allocator<_Ty>::_Not_user_specialized>
		: true_type
	{	
	};

		
struct _Alloc_allocate
	{	
		

	template<class _Alloc,
		class _Size_type,
		class _Const_void_pointer>
		static auto _Fn(int, _Alloc& _Al,
			_Size_type _Count,
			_Const_void_pointer _Hint)
			-> decltype(_Al.allocate(_Count, _Hint))
		{	
		return (_Al.allocate(_Count, _Hint));
		}

	template<class _Alloc,
		class _Size_type,
		class _Const_void_pointer>
		static auto _Fn(_Wrap_int, _Alloc& _Al,
			_Size_type _Count,
			_Const_void_pointer)
			-> decltype(_Al.allocate(_Count))
		{	
		return (_Al.allocate(_Count));
		}
	};

		
struct _Has_no_alloc_construct_tag
	{	
	};

template<class _Void,
	class... _Types>
	struct _Has_no_alloc_construct
		: true_type
	{	
	};

template<class _Alloc,
	class _Ptr,
	class... _Args>
	struct _Has_no_alloc_construct<
		void_t<
			_Has_no_alloc_construct_tag,	
			decltype(::std:: declval<_Alloc&>().construct(::std:: declval<_Ptr>(), ::std:: declval<_Args>()...))>,
		_Alloc, _Ptr, _Args...>
		: false_type
	{	
	};

template<class _Alloc,
	class _Ptr,
	class... _Args>
	using _Uses_default_construct = disjunction<
		_Is_default_allocator<_Alloc>,
		_Has_no_alloc_construct<void, _Alloc, _Ptr, _Args...>>;

template<class _Alloc,
	class _Ptr,
	class... _Args>
	using _Uses_default_construct_t = typename _Uses_default_construct<_Alloc, _Ptr, _Args...>::type;


		
struct _Has_no_alloc_destroy_tag
	{	
	};

template<class _Alloc,
	class _Ptr,
	class = void>
	struct _Has_no_alloc_destroy
		: true_type
	{	
	};

template<class _Alloc,
	class _Ptr>
	struct _Has_no_alloc_destroy<_Alloc, _Ptr, void_t<
			_Has_no_alloc_destroy_tag,	
			decltype(::std:: declval<_Alloc&>().destroy(::std:: declval<_Ptr>()))>>
		: false_type
	{	
	};

template<class _Alloc,
	class _Ptr>
	using _Uses_default_destroy = disjunction<
		_Is_default_allocator<_Alloc>,
		_Has_no_alloc_destroy<_Alloc, _Ptr>>;

template<class _Alloc,
	class _Ptr>
	using _Uses_default_destroy_t = typename _Uses_default_destroy<_Alloc, _Ptr>::type;


		
struct _Alloc_max_size
	{	
	template<class _Ty>
		static auto _Fn(int, const _Ty& _Al) noexcept
			-> decltype(_Al.max_size())
		{	
		return (_Al.max_size());
		}

	template<class _Ty>
		static auto _Fn(_Wrap_int, const _Ty&) noexcept
			-> typename _Get_size_type<_Ty>::type
		{	
		return ((numeric_limits<typename _Get_size_type<_Ty>::type>::max)()
			/ sizeof(typename _Ty::value_type));
		}
	};

		
struct _Alloc_select
	{	
		

	template<class _Ty>
		static auto _Fn(int, const _Ty& _Al)
			-> decltype((_Ty)_Al.select_on_container_copy_construction())
		{	
		return (_Al.select_on_container_copy_construction());
		}

	template<class _Ty>
		static auto _Fn(_Wrap_int, const _Ty& _Al)
			-> _Ty
		{	
		return (_Al);
		}
	};

		
template<class _Alloc>
	struct allocator_traits
	{	
	typedef _Alloc allocator_type;
	typedef typename _Alloc::value_type value_type;

	typedef typename _Get_pointer_type<_Alloc>::type
		pointer;
	typedef typename _Get_const_pointer_type<_Alloc>::type
		const_pointer;
	typedef typename _Get_void_pointer_type<_Alloc>::type
		void_pointer;
	typedef typename _Get_const_void_pointer_type<_Alloc>::type
		const_void_pointer;

	typedef typename _Get_size_type<_Alloc>::type size_type;
	typedef typename _Get_difference_type<_Alloc>::type difference_type;

	typedef typename _Get_propagate_on_container_copy<_Alloc>::type
		propagate_on_container_copy_assignment;
	typedef typename _Get_propagate_on_container_move<_Alloc>::type
		propagate_on_container_move_assignment;
	typedef typename _Get_propagate_on_container_swap<_Alloc>::type
		propagate_on_container_swap;
	typedef typename _Get_is_always_equal<_Alloc>::type
		is_always_equal;

	template<class _Other>
		using rebind_alloc = typename _Get_rebind_type<_Alloc, _Other>::type;

	template<class _Other>
		using rebind_traits = allocator_traits<rebind_alloc<_Other> >;

	static __declspec(allocator) pointer allocate(_Alloc& _Al, size_type _Count)
		{	
		return (_Al.allocate(_Count));
		}

	static __declspec(allocator) pointer allocate(_Alloc& _Al, size_type _Count,
		const_void_pointer _Hint)
		{	
		return (_Alloc_allocate::_Fn(0, _Al, _Count, _Hint));
		}

	static void deallocate(_Alloc& _Al,
		pointer _Ptr, size_type _Count)
		{	
		_Al.deallocate(_Ptr, _Count);
		}

	template<class _Ty,
		class... _Types>
		static void _Construct1(true_type, _Alloc&, _Ty *_Ptr,
			_Types&&... _Args)
		{	
		::new (static_cast<void *>(_Ptr))
			_Ty(::std:: forward<_Types>(_Args)...);
		}

	template<class _Ty,
		class... _Types>
		static void _Construct1(false_type, _Alloc& _Al, _Ty *_Ptr,
			_Types&&... _Args)
		{	
		_Al.construct(_Ptr, ::std:: forward<_Types>(_Args)...);
		}

	template<class _Ty,
		class... _Types>
		static void construct(_Alloc& _Al, _Ty *_Ptr,
			_Types&&... _Args)
		{	
		_Construct1(_Uses_default_construct_t<_Unwrap_alloc_t<_Alloc>, _Ty *, _Types...>(),
			_Al, _Ptr, ::std:: forward<_Types>(_Args)...);
		}

	template<class _Ty>
		static void _Destroy1(_Alloc&, _Ty *_Ptr, true_type)
		{	
		_Ptr->~_Ty();
		}

	template<class _Ty>
		static void _Destroy1(_Alloc& _Al, _Ty *_Ptr, false_type)
		{	
		_Al.destroy(_Ptr);
		}

	template<class _Ty>
		static void destroy(_Alloc& _Al, _Ty *_Ptr)
		{	
		_Destroy1(_Al, _Ptr, _Uses_default_destroy_t<_Unwrap_alloc_t<_Alloc>, _Ty *>());
		}

	static size_type max_size(const _Alloc& _Al) noexcept
		{	
		return (_Alloc_max_size::_Fn(0, _Al));
		}

	static _Alloc select_on_container_copy_construction(
		const _Alloc& _Al)
		{	
		return (_Alloc_select::_Fn(0, _Al));
		}
	};

		
template<class _Ty>
	class allocator
	{	
public:
	static_assert(!is_const<_Ty>::value,
		"The C++ Standard forbids containers of const elements "
		"because allocator<const T> is ill-formed.");

	typedef void _Not_user_specialized;

	typedef _Ty value_type;

	typedef value_type *pointer;
	typedef const value_type *const_pointer;

	typedef value_type& reference;
	typedef const value_type& const_reference;

	typedef size_t size_type;
	typedef ptrdiff_t difference_type;

	typedef true_type propagate_on_container_move_assignment;
	typedef true_type is_always_equal;

	template<class _Other>
		struct rebind
		{	
		typedef allocator<_Other> other;
		};

	pointer address(reference _Val) const noexcept
		{	
		return (::std:: addressof(_Val));
		}

	const_pointer address(const_reference _Val) const noexcept
		{	
		return (::std:: addressof(_Val));
		}

	allocator() noexcept
		{	
		}

	allocator(const allocator<_Ty>&) noexcept
		{	
		}

	template<class _Other>
		allocator(const allocator<_Other>&) noexcept
		{	
		}

	template<class _Other>
		allocator<_Ty>& operator=(const allocator<_Other>&)
		{	
		return (*this);
		}

	void deallocate(pointer _Ptr, size_type _Count)
		{	
		_Deallocate(_Ptr, _Count, sizeof (_Ty));
		}

	__declspec(allocator) pointer allocate(size_type _Count)
		{	
		return (static_cast<pointer>(_Allocate(_Count, sizeof (_Ty))));
		}

	__declspec(allocator) pointer allocate(size_type _Count, const void *)
		{	
		return (allocate(_Count));
		}

	template<class _Objty,
		class... _Types>
		void construct(_Objty *_Ptr, _Types&&... _Args)
		{	
		::new ((void *)_Ptr) _Objty(::std:: forward<_Types>(_Args)...);
		}


	template<class _Uty>
		void destroy(_Uty *_Ptr)
		{	
		_Ptr->~_Uty();
		}

	size_t max_size() const noexcept
		{	
		return ((size_t)(-1) / sizeof (_Ty));
		}
	};

		
template<>
	class allocator<void>
	{	
public:
	typedef void _Not_user_specialized;

	typedef void value_type;

	typedef void *pointer;
	typedef const void *const_pointer;

	template<class _Other>
		struct rebind
		{	
		typedef allocator<_Other> other;
		};

	allocator() noexcept
		{	
		}

	allocator(const allocator<void>&) noexcept
		{	
		}

	template<class _Other>
		allocator(const allocator<_Other>&) noexcept
		{	
		}

	template<class _Other>
		allocator<void>& operator=(const allocator<_Other>&)
		{	
		return (*this);
		}
	};

template<class _Ty,
	class _Other> inline
	bool operator==(const allocator<_Ty>&,
		const allocator<_Other>&) noexcept
	{	
	return (true);
	}

template<class _Ty,
	class _Other> inline
	bool operator!=(const allocator<_Ty>& _Left,
		const allocator<_Other>& _Right) noexcept
	{	
	return (false);
	}

		
template<class _Ty>
	struct allocator_traits<allocator<_Ty> >
	{	
	typedef allocator<_Ty> _Alloc;

	typedef _Alloc allocator_type;
	typedef _Ty value_type;

	typedef value_type *pointer;
	typedef const value_type *const_pointer;
	typedef void *void_pointer;
	typedef const void *const_void_pointer;

	typedef size_t size_type;
	typedef ptrdiff_t difference_type;

	typedef false_type propagate_on_container_copy_assignment;
	typedef true_type propagate_on_container_move_assignment;
	typedef false_type propagate_on_container_swap;
	typedef true_type is_always_equal;

	template<class _Other>
		using rebind_alloc = allocator<_Other>;

	template<class _Other>
		using rebind_traits = allocator_traits<allocator<_Other> >;

	static __declspec(allocator) pointer allocate(_Alloc& _Al, size_type _Count)
		{	
		return (_Al.allocate(_Count));
		}

	static __declspec(allocator) pointer allocate(_Alloc& _Al, size_type _Count,
		const_void_pointer _Hint)
		{	
		return (_Al.allocate(_Count, _Hint));
		}

	static void deallocate(_Alloc& _Al,
		pointer _Ptr, size_type _Count)
		{	
		_Al.deallocate(_Ptr, _Count);
		}

	template<class _Objty,
		class... _Types>
		static void construct(_Alloc& _Al, _Objty *_Ptr,
			_Types&&... _Args)
		{	
		_Al.construct(_Ptr, ::std:: forward<_Types>(_Args)...);
		}


	template<class _Uty>
		static void destroy(_Alloc& _Al, _Uty *_Ptr)
		{	
		_Al.destroy(_Ptr);
		}

	static size_type max_size(const _Alloc& _Al) noexcept
		{	
		return (_Al.max_size());
		}

	static _Alloc select_on_container_copy_construction(
		const _Alloc& _Al)
		{	
		return (_Al);
		}
	};

		
template<class _Alloc>
	struct _Wrap_alloc
		: public _Alloc
	{	
	typedef _Alloc _Mybase;
	typedef allocator_traits<_Alloc> _Mytraits;

	typedef typename _Mytraits::value_type value_type;

	typedef typename _Mytraits::pointer pointer;
	typedef typename _Mytraits::const_pointer const_pointer;
	typedef typename _Mytraits::void_pointer void_pointer;
	typedef typename _Mytraits::const_void_pointer const_void_pointer;

	typedef typename _If<is_void<value_type>::value,
		int, value_type>::type& reference;
	typedef typename _If<is_void<const value_type>::value,
		const int, const value_type>::type& const_reference;

	typedef typename _Mytraits::size_type size_type;
	typedef typename _Mytraits::difference_type difference_type;

	typedef typename _Mytraits::propagate_on_container_copy_assignment
		propagate_on_container_copy_assignment;
	typedef typename _Mytraits::propagate_on_container_move_assignment
		propagate_on_container_move_assignment;
	typedef typename _Mytraits::propagate_on_container_swap
		propagate_on_container_swap;
	typedef typename _Mytraits::is_always_equal
		is_always_equal;

	_Wrap_alloc select_on_container_copy_construction(_Nil = _Nil()) const
		{	
		return (_Mytraits::select_on_container_copy_construction(*this));
		}

	template<class _Other>
		struct rebind
		{	
		typedef typename _Mytraits::template rebind_alloc<_Other>
			_Other_alloc;
		typedef _Wrap_alloc<_Other_alloc> other;
		};

	pointer address(reference _Val) const
		{	
		return (pointer_traits<pointer>::pointer_to(_Val));
		}

	const_pointer address(const_reference _Val) const
		{	
		return (pointer_traits<const_pointer>::pointer_to(_Val));
		}

	_Wrap_alloc() noexcept(is_nothrow_default_constructible<_Alloc>::value)
		: _Mybase()
		{	
		}

	_Wrap_alloc(const _Wrap_alloc& _Right) noexcept
		: _Mybase(_Right)
		{	
		}

	_Wrap_alloc(_Wrap_alloc&& _Right) noexcept
		: _Mybase(::std:: move(_Right))
		{	
		}

	template<class _Other>
		_Wrap_alloc(_Other&& _Right) noexcept
		: _Mybase(::std:: forward<_Other>(_Right))
		{	
		}

	_Wrap_alloc& operator=(const _Wrap_alloc& _Right)
		{	
		_Mybase::operator=(_Right);
		return (*this);
		}

	_Wrap_alloc& operator=(_Wrap_alloc&& _Right)
		{	
		_Mybase::operator=(::std:: move(_Right));
		return (*this);
		}

	template<class _Other>
		_Wrap_alloc& operator=(_Other&& _Right)
		{	
		_Mybase::operator=(::std:: forward<_Other>(_Right));
		return (*this);
		}

	__declspec(allocator) pointer allocate(size_type _Count)
		{	
		return (_Mybase::allocate(_Count));
		}

	__declspec(allocator) pointer allocate(size_type _Count,
		const_void_pointer _Hint, _Nil = _Nil())
		{	
		return (_Mytraits::allocate(*this, _Count, _Hint));
		}

	void deallocate(pointer _Ptr, size_type _Count)
		{	
		_Mybase::deallocate(_Ptr, _Count);
		}

	template<class _Ty,
		class... _Types>
		void construct(_Ty *_Ptr,
			_Types&&... _Args)
		{	
		_Mytraits::construct(*this, _Ptr,
			::std:: forward<_Types>(_Args)...);
		}


	template<class _Ty>
		void destroy(_Ty *_Ptr)
		{	
		_Mytraits::destroy(*this, _Ptr);
		}

	size_type max_size(_Nil = _Nil()) const noexcept
		{	
		return (_Mytraits::max_size(*this));
		}
	};

template<class _Ty,
	class _Other> inline
	bool operator==(const _Wrap_alloc<_Ty>& _Left,
		const _Wrap_alloc<_Other>& _Right) noexcept
	{	
	return (static_cast<const _Ty&>(_Left)
		== static_cast<const _Other&>(_Right));
	}

template<class _Ty,
	class _Other> inline
	bool operator!=(const _Wrap_alloc<_Ty>& _Left,
		const _Wrap_alloc<_Other>& _Right) noexcept
	{	
	return (!(_Left == _Right));
	}

		
template<class _Alty> inline
	void _Pocca(_Alty& _Left, const _Alty& _Right, true_type) noexcept
	{	
	_Left = _Right;
	}

template<class _Alty> inline
	void _Pocca(_Alty&, const _Alty&, false_type) noexcept
	{	
	}

template<class _Alty> inline
	void _Pocca(_Alty& _Left, const _Alty& _Right) noexcept
	{	
	typename _Alty::propagate_on_container_copy_assignment _Tag;
	_Pocca(_Left, _Right, _Tag);
	}

		
template<class _Alty> inline
	void _Pocma(_Alty& _Left, _Alty& _Right, true_type) noexcept
	{	
	_Left = ::std:: move(_Right);
	}

template<class _Alty> inline
	void _Pocma(_Alty&, _Alty&, false_type) noexcept
	{	
	}

template<class _Alty> inline
	void _Pocma(_Alty& _Left, _Alty& _Right) noexcept
	{	
	typename _Alty::propagate_on_container_move_assignment _Tag;
	_Pocma(_Left, _Right, _Tag);
	}

		
template<class _Alty> inline
	void _Pocs(_Alty& _Left, _Alty& _Right, true_type) noexcept
	{	
	_Swap_adl(_Left, _Right);
	}

template<class _Alty> inline
	void _Pocs(_Alty& _Left, _Alty& _Right, false_type) noexcept
	{	
	if (_Left != _Right)
		{	
 


		::std:: terminate();
 
		}
	}

template<class _Alty> inline
	void _Pocs(_Alty& _Left, _Alty& _Right) noexcept
	{	
	typename _Alty::propagate_on_container_swap _Tag;
	_Pocs(_Left, _Right, _Tag);
	}


		
template<class _Alloc,
	class _Ptr = typename _Wrap_alloc<_Alloc>::pointer> inline
	void _Destroy_range1(_Ptr _First, _Ptr _Last, _Wrap_alloc<_Alloc>& _Al, false_type)
	{	
	for (; _First != _Last; ++_First)
		_Al.destroy(_Unfancy(_First));
	}

template<class _Alloc,
	class _Ptr = typename _Wrap_alloc<_Alloc>::pointer> inline
	void _Destroy_range1(_Ptr, _Ptr, _Wrap_alloc<_Alloc>&, true_type)
	{	
		
	}

template<class _Alloc,
	class _Ptr = typename _Wrap_alloc<_Alloc>::pointer> inline
	void _Destroy_range(_Ptr _First, _Ptr _Last, _Wrap_alloc<_Alloc>& _Al)
	{	
		
		
	typedef typename _Alloc::value_type _Val;
	_Destroy_range1(_First, _Last, _Al, typename conjunction<
		is_trivially_destructible<_Val>,
		_Uses_default_destroy<_Alloc, _Val *>>::type());
	}


		
template<class _FwdIt> inline
	void _Destroy_range1(_FwdIt _First, _FwdIt _Last, false_type)
	{	
	for (; _First != _Last; ++_First)
		_Destroy(_First);
	}

template<class _FwdIt> inline
	void _Destroy_range1(_FwdIt, _FwdIt, true_type)
	{	
		
	}

template<class _FwdIt> inline
	void _Destroy_range(_FwdIt _First, _FwdIt _Last)
	{	
		
		
	_Destroy_range1(_First, _Last, is_trivially_destructible<_Iter_value_t<_FwdIt>>());
	}
}

		
  

#pragma once





 #pragma pack(push,8)
 #pragma warning(push,3)
 
 

namespace std {
		
typedef enum memory_order {
	memory_order_relaxed,
	memory_order_consume,
	memory_order_acquire,
	memory_order_release,
	memory_order_acq_rel,
	memory_order_seq_cst
	} memory_order;

typedef _Uint32t _Uint4_t;
typedef _Uint4_t _Atomic_integral_t;

	
	




  
  
  
  
  

  
   
  



		

typedef long _Atomic_flag_t;

  

		
typedef _Atomic_integral_t _Atomic_counter_t;

inline _Atomic_integral_t
	_Get_atomic_count(const _Atomic_counter_t& _Counter)
	{	
	return (_Counter);
	}

inline void _Init_atomic_counter(_Atomic_counter_t& _Counter,
	_Atomic_integral_t _Value)
	{	
	_Counter = _Value;
	}

 
  
   
  


 

		
extern "C" {
__declspec(dllimport) void __cdecl _Lock_shared_ptr_spin_lock();
__declspec(dllimport) void __cdecl _Unlock_shared_ptr_spin_lock();
}
}
 
 #pragma warning(pop)
 #pragma pack(pop)









  
   












#pragma once













































































































































































































































































































































#pragma once












































































































































































































































































































































__pragma(pack(push, 8)) extern "C" {


























    typedef struct __declspec(align(16)) _SETJMP_FLOAT128
    {
        unsigned __int64 Part[2];
    } SETJMP_FLOAT128;

    
    typedef SETJMP_FLOAT128 _JBTYPE;

    typedef struct _JUMP_BUFFER
    {
        unsigned __int64 Frame;
        unsigned __int64 Rbx;
        unsigned __int64 Rsp;
        unsigned __int64 Rbp;
        unsigned __int64 Rsi;
        unsigned __int64 Rdi;
        unsigned __int64 R12;
        unsigned __int64 R13;
        unsigned __int64 R14;
        unsigned __int64 R15;
        unsigned __int64 Rip;
        unsigned long MxCsr;
        unsigned short FpCsr;
        unsigned short Spare;

        SETJMP_FLOAT128 Xmm6;
        SETJMP_FLOAT128 Xmm7;
        SETJMP_FLOAT128 Xmm8;
        SETJMP_FLOAT128 Xmm9;
        SETJMP_FLOAT128 Xmm10;
        SETJMP_FLOAT128 Xmm11;
        SETJMP_FLOAT128 Xmm12;
        SETJMP_FLOAT128 Xmm13;
        SETJMP_FLOAT128 Xmm14;
        SETJMP_FLOAT128 Xmm15;
    } _JUMP_BUFFER;




























































    
    typedef _JBTYPE jmp_buf[16];





    





int __cdecl _setjmp(
      jmp_buf _Buf
    );


    #pragma warning(push)
    #pragma warning(disable:4987) 
    __declspec(noreturn) void __cdecl longjmp(
          jmp_buf _Buf,
          int     _Value
        ) throw(...);
    #pragma warning(pop)








} __pragma(pack(pop))




    
        













#pragma once






























#pragma once































#pragma once































#pragma once
























#pragma once






























#pragma once





































#pragma once
















































#pragma once
































#pragma once













extern "C" { 




typedef union __declspec(intrin_type) __declspec(align(8)) __m64
{
    unsigned __int64    m64_u64;
    float               m64_f32[2];
    __int8              m64_i8[8];
    __int16             m64_i16[4];
    __int32             m64_i32[2];
    __int64             m64_i64;
    unsigned __int8     m64_u8[8];
    unsigned __int16    m64_u16[4];
    unsigned __int32    m64_u32[2];
} __m64;












































































































































}; 
























typedef union __declspec(intrin_type) __declspec(align(16)) __m128 {
     float               m128_f32[4];
     unsigned __int64    m128_u64[2];
     __int8              m128_i8[16];
     __int16             m128_i16[8];
     __int32             m128_i32[4];
     __int64             m128_i64[2];
     unsigned __int8     m128_u8[16];
     unsigned __int16    m128_u16[8];
     unsigned __int32    m128_u32[4];
 } __m128;







 
 
 
 
 
 
 
 
 
 




 
 
 
 
 
 
 
 
 
 












































































 
 
 


extern "C" { 
  






extern __m128 _mm_add_ss(__m128 _A, __m128 _B);
extern __m128 _mm_add_ps(__m128 _A, __m128 _B);
extern __m128 _mm_sub_ss(__m128 _A, __m128 _B);
extern __m128 _mm_sub_ps(__m128 _A, __m128 _B);
extern __m128 _mm_mul_ss(__m128 _A, __m128 _B);
extern __m128 _mm_mul_ps(__m128 _A, __m128 _B);
extern __m128 _mm_div_ss(__m128 _A, __m128 _B);
extern __m128 _mm_div_ps(__m128 _A, __m128 _B);
extern __m128 _mm_sqrt_ss(__m128 _A);
extern __m128 _mm_sqrt_ps(__m128 _A);
extern __m128 _mm_rcp_ss(__m128 _A);
extern __m128 _mm_rcp_ps(__m128 _A);
extern __m128 _mm_rsqrt_ss(__m128 _A);
extern __m128 _mm_rsqrt_ps(__m128 _A);
extern __m128 _mm_min_ss(__m128 _A, __m128 _B);
extern __m128 _mm_min_ps(__m128 _A, __m128 _B);
extern __m128 _mm_max_ss(__m128 _A, __m128 _B);
extern __m128 _mm_max_ps(__m128 _A, __m128 _B);





extern __m128 _mm_and_ps(__m128 _A, __m128 _B);
extern __m128 _mm_andnot_ps(__m128 _A, __m128 _B);
extern __m128 _mm_or_ps(__m128 _A, __m128 _B);
extern __m128 _mm_xor_ps(__m128 _A, __m128 _B);





extern __m128 _mm_cmpeq_ss(__m128 _A, __m128 _B);
extern __m128 _mm_cmpeq_ps(__m128 _A, __m128 _B);
extern __m128 _mm_cmplt_ss(__m128 _A, __m128 _B);
extern __m128 _mm_cmplt_ps(__m128 _A, __m128 _B);
extern __m128 _mm_cmple_ss(__m128 _A, __m128 _B);
extern __m128 _mm_cmple_ps(__m128 _A, __m128 _B);
extern __m128 _mm_cmpgt_ss(__m128 _A, __m128 _B);
extern __m128 _mm_cmpgt_ps(__m128 _A, __m128 _B);
extern __m128 _mm_cmpge_ss(__m128 _A, __m128 _B);
extern __m128 _mm_cmpge_ps(__m128 _A, __m128 _B);
extern __m128 _mm_cmpneq_ss(__m128 _A, __m128 _B);
extern __m128 _mm_cmpneq_ps(__m128 _A, __m128 _B);
extern __m128 _mm_cmpnlt_ss(__m128 _A, __m128 _B);
extern __m128 _mm_cmpnlt_ps(__m128 _A, __m128 _B);
extern __m128 _mm_cmpnle_ss(__m128 _A, __m128 _B);
extern __m128 _mm_cmpnle_ps(__m128 _A, __m128 _B);
extern __m128 _mm_cmpngt_ss(__m128 _A, __m128 _B);
extern __m128 _mm_cmpngt_ps(__m128 _A, __m128 _B);
extern __m128 _mm_cmpnge_ss(__m128 _A, __m128 _B);
extern __m128 _mm_cmpnge_ps(__m128 _A, __m128 _B);
extern __m128 _mm_cmpord_ss(__m128 _A, __m128 _B);
extern __m128 _mm_cmpord_ps(__m128 _A, __m128 _B);
extern __m128 _mm_cmpunord_ss(__m128 _A, __m128 _B);
extern __m128 _mm_cmpunord_ps(__m128 _A, __m128 _B);
extern int _mm_comieq_ss(__m128 _A, __m128 _B);
extern int _mm_comilt_ss(__m128 _A, __m128 _B);
extern int _mm_comile_ss(__m128 _A, __m128 _B);
extern int _mm_comigt_ss(__m128 _A, __m128 _B);
extern int _mm_comige_ss(__m128 _A, __m128 _B);
extern int _mm_comineq_ss(__m128 _A, __m128 _B);
extern int _mm_ucomieq_ss(__m128 _A, __m128 _B);
extern int _mm_ucomilt_ss(__m128 _A, __m128 _B);
extern int _mm_ucomile_ss(__m128 _A, __m128 _B);
extern int _mm_ucomigt_ss(__m128 _A, __m128 _B);
extern int _mm_ucomige_ss(__m128 _A, __m128 _B);
extern int _mm_ucomineq_ss(__m128 _A, __m128 _B);





extern int _mm_cvt_ss2si(__m128 _A);
extern int _mm_cvtt_ss2si(__m128 _A);
extern __m128 _mm_cvt_si2ss(__m128, int);
extern float _mm_cvtss_f32(__m128 _A);














extern __int64 _mm_cvtss_si64(__m128 _A);
extern __int64 _mm_cvttss_si64(__m128 _A);
extern __m128  _mm_cvtsi64_ss(__m128 _A, __int64 _B);






extern __m128 _mm_shuffle_ps(__m128 _A, __m128 _B, unsigned int _Imm8);
extern __m128 _mm_unpackhi_ps(__m128 _A, __m128 _B);
extern __m128 _mm_unpacklo_ps(__m128 _A, __m128 _B);
extern __m128 _mm_loadh_pi(__m128, __m64 const*);
extern __m128 _mm_movehl_ps(__m128, __m128);
extern __m128 _mm_movelh_ps(__m128, __m128);
extern void _mm_storeh_pi(__m64 *, __m128);
extern __m128 _mm_loadl_pi(__m128, __m64 const*);
extern void _mm_storel_pi(__m64 *, __m128);
extern int _mm_movemask_ps(__m128 _A);

























extern __m128 _mm_set_ss(float _A);
extern __m128 _mm_set_ps1(float _A);
extern __m128 _mm_set_ps(float _A, float _B, float _C, float _D);
extern __m128 _mm_setr_ps(float _A, float _B, float _C, float _D);
extern __m128 _mm_setzero_ps(void);
extern __m128 _mm_load_ss(float const*_A);
extern __m128 _mm_load_ps1(float const*_A);
extern __m128 _mm_load_ps(float const*_A);
extern __m128 _mm_loadr_ps(float const*_A);
extern __m128 _mm_loadu_ps(float const*_A);
extern void _mm_store_ss(float *_V, __m128 _A);
extern void _mm_store_ps1(float *_V, __m128 _A);
extern void _mm_store_ps(float *_V, __m128 _A);
extern void _mm_storer_ps(float *_V, __m128 _A);
extern void _mm_storeu_ps(float *_V, __m128 _A);
extern void _mm_prefetch(char const*_A, int _Sel);



extern void _mm_stream_ps(float *, __m128);
extern __m128 _mm_move_ss(__m128 _A, __m128 _B);

extern void _mm_sfence(void);
extern unsigned int _mm_getcsr(void);
extern void _mm_setcsr(unsigned int);
































 
 
 






















































































































}; 







typedef union __declspec(intrin_type) __declspec(align(16)) __m128i {
    __int8              m128i_i8[16];
    __int16             m128i_i16[8];
    __int32             m128i_i32[4];
    __int64             m128i_i64[2];
    unsigned __int8     m128i_u8[16];
    unsigned __int16    m128i_u16[8];
    unsigned __int32    m128i_u32[4];
    unsigned __int64    m128i_u64[2];
} __m128i;

typedef struct __declspec(intrin_type) __declspec(align(16)) __m128d {
    double              m128d_f64[2];
} __m128d;






 
 
 


extern "C" { 
  






extern __m128d _mm_add_sd(__m128d _A, __m128d _B);
extern __m128d _mm_add_pd(__m128d _A, __m128d _B);
extern __m128d _mm_sub_sd(__m128d _A, __m128d _B);
extern __m128d _mm_sub_pd(__m128d _A, __m128d _B);
extern __m128d _mm_mul_sd(__m128d _A, __m128d _B);
extern __m128d _mm_mul_pd(__m128d _A, __m128d _B);
extern __m128d _mm_sqrt_sd(__m128d _A, __m128d _B);
extern __m128d _mm_sqrt_pd(__m128d _A);
extern __m128d _mm_div_sd(__m128d _A, __m128d _B);
extern __m128d _mm_div_pd(__m128d _A, __m128d _B);
extern __m128d _mm_min_sd(__m128d _A, __m128d _B);
extern __m128d _mm_min_pd(__m128d _A, __m128d _B);
extern __m128d _mm_max_sd(__m128d _A, __m128d _B);
extern __m128d _mm_max_pd(__m128d _A, __m128d _B);





extern __m128d _mm_and_pd(__m128d _A, __m128d _B);
extern __m128d _mm_andnot_pd(__m128d _A, __m128d _B);
extern __m128d _mm_or_pd(__m128d _A, __m128d _B);
extern __m128d _mm_xor_pd(__m128d _A, __m128d _B);





extern __m128d _mm_cmpeq_sd(__m128d _A, __m128d _B);
extern __m128d _mm_cmpeq_pd(__m128d _A, __m128d _B);
extern __m128d _mm_cmplt_sd(__m128d _A, __m128d _B);
extern __m128d _mm_cmplt_pd(__m128d _A, __m128d _B);
extern __m128d _mm_cmple_sd(__m128d _A, __m128d _B);
extern __m128d _mm_cmple_pd(__m128d _A, __m128d _B);
extern __m128d _mm_cmpgt_sd(__m128d _A, __m128d _B);
extern __m128d _mm_cmpgt_pd(__m128d _A, __m128d _B);
extern __m128d _mm_cmpge_sd(__m128d _A, __m128d _B);
extern __m128d _mm_cmpge_pd(__m128d _A, __m128d _B);
extern __m128d _mm_cmpneq_sd(__m128d _A, __m128d _B);
extern __m128d _mm_cmpneq_pd(__m128d _A, __m128d _B);
extern __m128d _mm_cmpnlt_sd(__m128d _A, __m128d _B);
extern __m128d _mm_cmpnlt_pd(__m128d _A, __m128d _B);
extern __m128d _mm_cmpnle_sd(__m128d _A, __m128d _B);
extern __m128d _mm_cmpnle_pd(__m128d _A, __m128d _B);
extern __m128d _mm_cmpngt_sd(__m128d _A, __m128d _B);
extern __m128d _mm_cmpngt_pd(__m128d _A, __m128d _B);
extern __m128d _mm_cmpnge_sd(__m128d _A, __m128d _B);
extern __m128d _mm_cmpnge_pd(__m128d _A, __m128d _B);
extern __m128d _mm_cmpord_pd(__m128d _A, __m128d _B);
extern __m128d _mm_cmpord_sd(__m128d _A, __m128d _B);
extern __m128d _mm_cmpunord_pd(__m128d _A, __m128d _B);
extern __m128d _mm_cmpunord_sd(__m128d _A, __m128d _B);
extern int _mm_comieq_sd(__m128d _A, __m128d _B);
extern int _mm_comilt_sd(__m128d _A, __m128d _B);
extern int _mm_comile_sd(__m128d _A, __m128d _B);
extern int _mm_comigt_sd(__m128d _A, __m128d _B);
extern int _mm_comige_sd(__m128d _A, __m128d _B);
extern int _mm_comineq_sd(__m128d _A, __m128d _B);
extern int _mm_ucomieq_sd(__m128d _A, __m128d _B);
extern int _mm_ucomilt_sd(__m128d _A, __m128d _B);
extern int _mm_ucomile_sd(__m128d _A, __m128d _B);
extern int _mm_ucomigt_sd(__m128d _A, __m128d _B);
extern int _mm_ucomige_sd(__m128d _A, __m128d _B);
extern int _mm_ucomineq_sd(__m128d _A, __m128d _B);





extern __m128d _mm_cvtepi32_pd(__m128i _A);
extern __m128i _mm_cvtpd_epi32(__m128d _A);
extern __m128i _mm_cvttpd_epi32(__m128d _A);
extern __m128 _mm_cvtepi32_ps(__m128i _A);
extern __m128i _mm_cvtps_epi32(__m128 _A);
extern __m128i _mm_cvttps_epi32(__m128 _A);
extern __m128 _mm_cvtpd_ps(__m128d _A);
extern __m128d _mm_cvtps_pd(__m128 _A);
extern __m128 _mm_cvtsd_ss(__m128 _A, __m128d _B);
extern __m128d _mm_cvtss_sd(__m128d _A, __m128 _B);

extern int _mm_cvtsd_si32(__m128d _A);
extern int _mm_cvttsd_si32(__m128d _A);
extern __m128d _mm_cvtsi32_sd(__m128d _A, int _B);











extern __m128d _mm_unpackhi_pd(__m128d _A, __m128d _B);
extern __m128d _mm_unpacklo_pd(__m128d _A, __m128d _B);
extern int _mm_movemask_pd(__m128d _A);
extern __m128d _mm_shuffle_pd(__m128d _A, __m128d _B, int _I);





extern __m128d _mm_load_pd(double const*_Dp);
extern __m128d _mm_load1_pd(double const*_Dp);
extern __m128d _mm_loadr_pd(double const*_Dp);
extern __m128d _mm_loadu_pd(double const*_Dp);
extern __m128d _mm_load_sd(double const*_Dp);
extern __m128d _mm_loadh_pd(__m128d _A, double const*_Dp);
extern __m128d _mm_loadl_pd(__m128d _A, double const*_Dp);





extern __m128d _mm_set_sd(double _W);
extern __m128d _mm_set1_pd(double _A);
extern __m128d _mm_set_pd(double _Z, double _Y);
extern __m128d _mm_setr_pd(double _Y, double _Z);
extern __m128d _mm_setzero_pd(void);
extern __m128d _mm_move_sd(__m128d _A, __m128d _B);





extern void _mm_store_sd(double *_Dp, __m128d _A);
extern void _mm_store1_pd(double *_Dp, __m128d _A);
extern void _mm_store_pd(double *_Dp, __m128d _A);
extern void _mm_storeu_pd(double *_Dp, __m128d _A);
extern void _mm_storer_pd(double *_Dp, __m128d _A);
extern void _mm_storeh_pd(double *_Dp, __m128d _A);
extern void _mm_storel_pd(double *_Dp, __m128d _A);





extern __m128i _mm_add_epi8(__m128i _A, __m128i _B);
extern __m128i _mm_add_epi16(__m128i _A, __m128i _B);
extern __m128i _mm_add_epi32(__m128i _A, __m128i _B);



extern __m128i _mm_add_epi64(__m128i _A, __m128i _B);
extern __m128i _mm_adds_epi8(__m128i _A, __m128i _B);
extern __m128i _mm_adds_epi16(__m128i _A, __m128i _B);
extern __m128i _mm_adds_epu8(__m128i _A, __m128i _B);
extern __m128i _mm_adds_epu16(__m128i _A, __m128i _B);
extern __m128i _mm_avg_epu8(__m128i _A, __m128i _B);
extern __m128i _mm_avg_epu16(__m128i _A, __m128i _B);
extern __m128i _mm_madd_epi16(__m128i _A, __m128i _B);
extern __m128i _mm_max_epi16(__m128i _A, __m128i _B);
extern __m128i _mm_max_epu8(__m128i _A, __m128i _B);
extern __m128i _mm_min_epi16(__m128i _A, __m128i _B);
extern __m128i _mm_min_epu8(__m128i _A, __m128i _B);
extern __m128i _mm_mulhi_epi16(__m128i _A, __m128i _B);
extern __m128i _mm_mulhi_epu16(__m128i _A, __m128i _B);
extern __m128i _mm_mullo_epi16(__m128i _A, __m128i _B);



extern __m128i _mm_mul_epu32(__m128i _A, __m128i _B);
extern __m128i _mm_sad_epu8(__m128i _A, __m128i _B);
extern __m128i _mm_sub_epi8(__m128i _A, __m128i _B);
extern __m128i _mm_sub_epi16(__m128i _A, __m128i _B);
extern __m128i _mm_sub_epi32(__m128i _A, __m128i _B);



extern __m128i _mm_sub_epi64(__m128i _A, __m128i _B);
extern __m128i _mm_subs_epi8(__m128i _A, __m128i _B);
extern __m128i _mm_subs_epi16(__m128i _A, __m128i _B);
extern __m128i _mm_subs_epu8(__m128i _A, __m128i _B);
extern __m128i _mm_subs_epu16(__m128i _A, __m128i _B);





extern __m128i _mm_and_si128(__m128i _A, __m128i _B);
extern __m128i _mm_andnot_si128(__m128i _A, __m128i _B);
extern __m128i _mm_or_si128(__m128i _A, __m128i _B);
extern __m128i _mm_xor_si128(__m128i _A, __m128i _B);





extern __m128i _mm_slli_si128(__m128i _A, int _Imm);
extern __m128i _mm_slli_epi16(__m128i _A, int _Count);
extern __m128i _mm_sll_epi16(__m128i _A, __m128i _Count);
extern __m128i _mm_slli_epi32(__m128i _A, int _Count);
extern __m128i _mm_sll_epi32(__m128i _A, __m128i _Count);
extern __m128i _mm_slli_epi64(__m128i _A, int _Count);
extern __m128i _mm_sll_epi64(__m128i _A, __m128i _Count);
extern __m128i _mm_srai_epi16(__m128i _A, int _Count);
extern __m128i _mm_sra_epi16(__m128i _A, __m128i _Count);
extern __m128i _mm_srai_epi32(__m128i _A, int _Count);
extern __m128i _mm_sra_epi32(__m128i _A, __m128i _Count);
extern __m128i _mm_srli_si128(__m128i _A, int _Imm);
extern __m128i _mm_srli_epi16(__m128i _A, int _Count);
extern __m128i _mm_srl_epi16(__m128i _A, __m128i _Count);
extern __m128i _mm_srli_epi32(__m128i _A, int _Count);
extern __m128i _mm_srl_epi32(__m128i _A, __m128i _Count);
extern __m128i _mm_srli_epi64(__m128i _A, int _Count);
extern __m128i _mm_srl_epi64(__m128i _A, __m128i _Count);





extern __m128i _mm_cmpeq_epi8(__m128i _A, __m128i _B);
extern __m128i _mm_cmpeq_epi16(__m128i _A, __m128i _B);
extern __m128i _mm_cmpeq_epi32(__m128i _A, __m128i _B);
extern __m128i _mm_cmpgt_epi8(__m128i _A, __m128i _B);
extern __m128i _mm_cmpgt_epi16(__m128i _A, __m128i _B);
extern __m128i _mm_cmpgt_epi32(__m128i _A, __m128i _B);
extern __m128i _mm_cmplt_epi8(__m128i _A, __m128i _B);
extern __m128i _mm_cmplt_epi16(__m128i _A, __m128i _B);
extern __m128i _mm_cmplt_epi32(__m128i _A, __m128i _B);





extern __m128i _mm_cvtsi32_si128(int _A);
extern int _mm_cvtsi128_si32(__m128i _A);





extern __m128i _mm_packs_epi16(__m128i _A, __m128i _B);
extern __m128i _mm_packs_epi32(__m128i _A, __m128i _B);
extern __m128i _mm_packus_epi16(__m128i _A, __m128i _B);
extern int _mm_extract_epi16(__m128i _A, int _Imm);
extern __m128i _mm_insert_epi16(__m128i _A, int _B, int _Imm);
extern int _mm_movemask_epi8(__m128i _A);
extern __m128i _mm_shuffle_epi32(__m128i _A, int _Imm);
extern __m128i _mm_shufflehi_epi16(__m128i _A, int _Imm);
extern __m128i _mm_shufflelo_epi16(__m128i _A, int _Imm);
extern __m128i _mm_unpackhi_epi8(__m128i _A, __m128i _B);
extern __m128i _mm_unpackhi_epi16(__m128i _A, __m128i _B);
extern __m128i _mm_unpackhi_epi32(__m128i _A, __m128i _B);
extern __m128i _mm_unpackhi_epi64(__m128i _A, __m128i _B);
extern __m128i _mm_unpacklo_epi8(__m128i _A, __m128i _B);
extern __m128i _mm_unpacklo_epi16(__m128i _A, __m128i _B);
extern __m128i _mm_unpacklo_epi32(__m128i _A, __m128i _B);
extern __m128i _mm_unpacklo_epi64(__m128i _A, __m128i _B);





extern __m128i _mm_load_si128(__m128i const*_P);
extern __m128i _mm_loadu_si128(__m128i const*_P);
extern __m128i _mm_loadl_epi64(__m128i const*_P);








extern __m128i _mm_set_epi64x(__int64 _I1,__int64 _I0);
extern __m128i _mm_set_epi32(int _I3, int _I2, int _I1, int _I0);
extern __m128i _mm_set_epi16(short _W7, short _W6, short _W5, short _W4,
                             short _W3, short _W2, short _W1, short _W0);
extern __m128i _mm_set_epi8(char _B15, char _B14, char _B13, char _B12,
                            char _B11, char _B10, char _B9, char _B8,
                            char _B7, char _B6, char _B5, char _B4,
                            char _B3, char _B2, char _B1, char _B0);



extern __m128i _mm_set1_epi64x(__int64 i);
extern __m128i _mm_set1_epi32(int _I);
extern __m128i _mm_set1_epi16(short _W);
extern __m128i _mm_set1_epi8(char _B);
extern __m128i _mm_setl_epi64(__m128i _Q);



extern __m128i _mm_setr_epi32(int _I0, int _I1, int _I2, int _I3);
extern __m128i _mm_setr_epi16(short _W0, short _W1, short _W2, short _W3,
                              short _W4, short _W5, short _W6, short _W7);
extern __m128i _mm_setr_epi8(char _B15, char _B14, char _B13, char _B12,
                             char _B11, char _B10, char _B9, char _B8,
                             char _B7, char _B6, char _B5, char _B4,
                             char _B3, char _B2, char _B1, char _B0);
extern __m128i _mm_setzero_si128(void);





extern void _mm_store_si128(__m128i *_P, __m128i _B);
extern void _mm_storeu_si128(__m128i *_P, __m128i _B);
extern void _mm_storel_epi64(__m128i *_P, __m128i _Q);
extern void _mm_maskmoveu_si128(__m128i _D, __m128i _N, char *_P);





extern __m128i _mm_move_epi64(__m128i _Q);









extern void _mm_stream_pd(double *_Dp, __m128d _A);
extern void _mm_stream_si128(__m128i *_P, __m128i _A);
extern void _mm_clflush(void const*_P);
extern void _mm_lfence(void);
extern void _mm_mfence(void);
extern void _mm_stream_si32(int *_P, int _I);
extern void _mm_pause(void);





extern double _mm_cvtsd_f64(__m128d _A);







extern __m128  _mm_castpd_ps(__m128d);
extern __m128i _mm_castpd_si128(__m128d);
extern __m128d _mm_castps_pd(__m128);
extern __m128i _mm_castps_si128(__m128);
extern __m128  _mm_castsi128_ps(__m128i);
extern __m128d _mm_castsi128_pd(__m128i);






extern __int64 _mm_cvtsd_si64(__m128d);
extern __int64 _mm_cvttsd_si64(__m128d);
extern __m128d _mm_cvtsi64_sd(__m128d, __int64);
extern __m128i _mm_cvtsi64_si128(__int64);
extern __int64 _mm_cvtsi128_si64(__m128i);





}; 







 
 
 














 
 
 


extern "C" { 
  






extern __m128 _mm_addsub_ps(__m128 , __m128 );
extern __m128 _mm_hadd_ps(__m128 , __m128 );
extern __m128 _mm_hsub_ps(__m128 , __m128 );
extern __m128 _mm_movehdup_ps(__m128 );
extern __m128 _mm_moveldup_ps(__m128 );





extern __m128d _mm_addsub_pd(__m128d , __m128d );
extern __m128d _mm_hadd_pd(__m128d , __m128d );
extern __m128d _mm_hsub_pd(__m128d , __m128d );
extern __m128d _mm_loaddup_pd(double const * );
extern __m128d _mm_movedup_pd(__m128d );




extern __m128i _mm_lddqu_si128(__m128i const * );







extern void _mm_monitor(void const * , unsigned , unsigned );




extern void _mm_mwait(unsigned , unsigned );


}; 















extern "C" {


    
    
    
    
    
    

    extern __m128i _mm_hadd_epi16 (__m128i, __m128i);
    extern __m128i _mm_hadd_epi32 (__m128i, __m128i);
    extern __m128i _mm_hadds_epi16 (__m128i, __m128i);







    
    
    
    
    
    
    

    extern __m128i _mm_hsub_epi16 (__m128i, __m128i);
    extern __m128i _mm_hsub_epi32 (__m128i, __m128i);
    extern __m128i _mm_hsubs_epi16 (__m128i, __m128i);







    
    
    
    
    
    
    
    

    extern __m128i _mm_maddubs_epi16 (__m128i, __m128i);





    
    

    extern __m128i _mm_mulhrs_epi16 (__m128i, __m128i);





    
    

    extern __m128i _mm_shuffle_epi8 (__m128i, __m128i);





    
    

    extern __m128i _mm_sign_epi8 (__m128i, __m128i);
    extern __m128i _mm_sign_epi16 (__m128i, __m128i);
    extern __m128i _mm_sign_epi32 (__m128i, __m128i);







    
    

    extern __m128i _mm_alignr_epi8 (__m128i, __m128i, int);





    
    

    extern __m128i _mm_abs_epi8 (__m128i);
    extern __m128i _mm_abs_epi16 (__m128i);
    extern __m128i _mm_abs_epi32 (__m128i);








};
























































extern "C" {


        
        

        extern __m128i _mm_blend_epi16 (__m128i, __m128i, const int );
        extern __m128i _mm_blendv_epi8 (__m128i, __m128i, __m128i mask);

        
        

        extern __m128  _mm_blend_ps (__m128, __m128, const int );
        extern __m128  _mm_blendv_ps(__m128, __m128, __m128 );

        
        

        extern __m128d _mm_blend_pd (__m128d, __m128d, const int );
        extern __m128d _mm_blendv_pd(__m128d, __m128d, __m128d );

        
        

        extern __m128  _mm_dp_ps(__m128, __m128, const int );
        extern __m128d _mm_dp_pd(__m128d, __m128d, const int );

        
        

        extern __m128i _mm_cmpeq_epi64(__m128i, __m128i);

        

        extern __m128i _mm_min_epi8 (__m128i, __m128i);
        extern __m128i _mm_max_epi8 (__m128i, __m128i);

        extern __m128i _mm_min_epu16(__m128i, __m128i);
        extern __m128i _mm_max_epu16(__m128i, __m128i);

        extern __m128i _mm_min_epi32(__m128i, __m128i);
        extern __m128i _mm_max_epi32(__m128i, __m128i);
        extern __m128i _mm_min_epu32(__m128i, __m128i);
        extern __m128i _mm_max_epu32(__m128i, __m128i);

        
        

        extern __m128i _mm_mullo_epi32(__m128i, __m128i);

        
        

        extern __m128i _mm_mul_epi32(__m128i, __m128i);

        
        

        extern int _mm_testz_si128(__m128i , __m128i );

        
        

        extern int _mm_testc_si128(__m128i , __m128i );

        
        
        

        extern int _mm_testnzc_si128(__m128i , __m128i );

        
        
        
        
        

        extern __m128 _mm_insert_ps(__m128 , __m128 , const int );

        




        
        

        extern int _mm_extract_ps(__m128 , const int );

        
        




        
        





        
        

        extern __m128i _mm_insert_epi8 (__m128i , int , const int );
        extern __m128i _mm_insert_epi32(__m128i , int , const int );


        extern __m128i _mm_insert_epi64(__m128i , __int64 , const int );

        
        

        extern int   _mm_extract_epi8 (__m128i , const int );
        extern int   _mm_extract_epi32(__m128i , const int );


        extern __int64 _mm_extract_epi64(__m128i , const int );


        
        

        extern __m128i _mm_minpos_epu16(__m128i);

        

        extern __m128d _mm_round_pd(__m128d , int );
        extern __m128d _mm_round_sd(__m128d , __m128d , int );

        

        extern __m128  _mm_round_ps(__m128  , int );
        extern __m128  _mm_round_ss(__m128 , __m128  , int );

        

        extern __m128i _mm_cvtepi8_epi32 (__m128i);
        extern __m128i _mm_cvtepi16_epi32(__m128i);
        extern __m128i _mm_cvtepi8_epi64 (__m128i);
        extern __m128i _mm_cvtepi32_epi64(__m128i);
        extern __m128i _mm_cvtepi16_epi64(__m128i);
        extern __m128i _mm_cvtepi8_epi16 (__m128i);

        

        extern __m128i _mm_cvtepu8_epi32 (__m128i);
        extern __m128i _mm_cvtepu16_epi32(__m128i);
        extern __m128i _mm_cvtepu8_epi64 (__m128i);
        extern __m128i _mm_cvtepu32_epi64(__m128i);
        extern __m128i _mm_cvtepu16_epi64(__m128i);
        extern __m128i _mm_cvtepu8_epi16 (__m128i);


        
        

        extern __m128i _mm_packus_epi32(__m128i, __m128i);

        
        
        

        extern __m128i _mm_mpsadbw_epu8(__m128i , __m128i , const int );

        



        extern __m128i _mm_stream_load_si128(const __m128i*);


}; 









extern "C" {














































    extern __m128i _mm_cmpistrm (__m128i , __m128i , const int );
    extern int     _mm_cmpistri (__m128i , __m128i , const int );

    extern __m128i _mm_cmpestrm (__m128i , int , __m128i , int , const int );
    extern int     _mm_cmpestri (__m128i , int , __m128i , int , const int );





    extern int     _mm_cmpistrz (__m128i , __m128i , const int );
    extern int     _mm_cmpistrc (__m128i , __m128i , const int );
    extern int     _mm_cmpistrs (__m128i , __m128i , const int );
    extern int     _mm_cmpistro (__m128i , __m128i , const int );
    extern int     _mm_cmpistra (__m128i , __m128i , const int );

    extern int     _mm_cmpestrz (__m128i , int , __m128i , int , const int );
    extern int     _mm_cmpestrc (__m128i , int , __m128i , int , const int );
    extern int     _mm_cmpestrs (__m128i , int , __m128i , int , const int );
    extern int     _mm_cmpestro (__m128i , int , __m128i , int , const int );
    extern int     _mm_cmpestra (__m128i , int , __m128i , int , const int );






    extern __m128i _mm_cmpgt_epi64(__m128i , __m128i );





    extern int _mm_popcnt_u32(unsigned int );


    extern __int64 _mm_popcnt_u64(unsigned __int64 );






    extern unsigned int _mm_crc32_u8 (unsigned int , unsigned char );
    extern unsigned int _mm_crc32_u16(unsigned int , unsigned short );
    extern unsigned int _mm_crc32_u32(unsigned int , unsigned int );


    extern unsigned __int64 _mm_crc32_u64(unsigned __int64 , unsigned __int64 );



}; 









extern "C" {






extern __m128i _mm_aesdec_si128(__m128i , __m128i );





extern __m128i _mm_aesdeclast_si128(__m128i , __m128i );





extern __m128i _mm_aesenc_si128(__m128i , __m128i );





extern __m128i _mm_aesenclast_si128(__m128i , __m128i );





extern __m128i _mm_aesimc_si128(__m128i );






extern __m128i _mm_aeskeygenassist_si128(__m128i , const int );







extern __m128i _mm_clmulepi64_si128(__m128i , __m128i ,
                                            const int );



}; 








extern "C" {





typedef union __declspec(intrin_type) __declspec(align(32)) __m256 {
    float m256_f32[8];
} __m256;

typedef struct __declspec(intrin_type) __declspec(align(32)) __m256d {
    double m256d_f64[4];
} __m256d;

typedef union  __declspec(intrin_type) __declspec(align(32)) __m256i {
    __int8              m256i_i8[32];
    __int16             m256i_i16[16];
    __int32             m256i_i32[8];
    __int64             m256i_i64[4];
    unsigned __int8     m256i_u8[32];
    unsigned __int16    m256i_u16[16];
    unsigned __int32    m256i_u32[8];
    unsigned __int64    m256i_u64[4];
} __m256i;



















































extern __m256d __cdecl _mm256_add_pd(__m256d, __m256d);









extern __m256 __cdecl _mm256_add_ps(__m256, __m256);












extern __m256d __cdecl _mm256_addsub_pd(__m256d, __m256d);












extern __m256 __cdecl _mm256_addsub_ps(__m256, __m256);








extern __m256d __cdecl _mm256_and_pd(__m256d, __m256d);








extern __m256 __cdecl _mm256_and_ps(__m256, __m256);








extern __m256d __cdecl _mm256_andnot_pd(__m256d, __m256d);








extern __m256 __cdecl _mm256_andnot_ps(__m256, __m256);













extern __m256d __cdecl _mm256_blend_pd(__m256d, __m256d, const int);













extern __m256 __cdecl _mm256_blend_ps(__m256, __m256, const int);









extern __m256d __cdecl _mm256_blendv_pd(__m256d, __m256d, __m256d);









extern __m256 __cdecl _mm256_blendv_ps(__m256, __m256, __m256);








extern __m256d __cdecl _mm256_div_pd(__m256d, __m256d);








extern __m256 __cdecl _mm256_div_ps(__m256, __m256);














extern __m256 __cdecl _mm256_dp_ps(__m256, __m256, const int);








extern __m256d __cdecl _mm256_hadd_pd(__m256d, __m256d);








extern __m256 __cdecl _mm256_hadd_ps(__m256, __m256);








extern __m256d __cdecl _mm256_hsub_pd(__m256d, __m256d);








extern __m256 __cdecl _mm256_hsub_ps(__m256, __m256);








extern __m256d __cdecl _mm256_max_pd(__m256d, __m256d);








extern __m256 __cdecl _mm256_max_ps(__m256, __m256);








extern __m256d __cdecl _mm256_min_pd(__m256d, __m256d);








extern __m256 __cdecl _mm256_min_ps(__m256, __m256);









extern __m256d __cdecl _mm256_mul_pd(__m256d, __m256d);









extern __m256 __cdecl _mm256_mul_ps(__m256, __m256);








extern __m256d __cdecl _mm256_or_pd(__m256d, __m256d);








extern __m256 __cdecl _mm256_or_ps(__m256, __m256);











extern __m256d __cdecl _mm256_shuffle_pd(__m256d, __m256d, const int);












extern __m256 __cdecl _mm256_shuffle_ps(__m256, __m256, const int);








extern __m256d __cdecl _mm256_sub_pd(__m256d, __m256d);









extern __m256 __cdecl _mm256_sub_ps(__m256, __m256);








extern __m256d __cdecl _mm256_xor_pd(__m256d, __m256d);








extern __m256 __cdecl _mm256_xor_ps(__m256, __m256);















extern __m128d __cdecl _mm_cmp_pd(__m128d, __m128d, const int);
extern __m256d __cdecl _mm256_cmp_pd(__m256d, __m256d, const int);















extern __m128 __cdecl _mm_cmp_ps(__m128, __m128, const int);
extern __m256 __cdecl _mm256_cmp_ps(__m256, __m256, const int);












extern __m128d __cdecl _mm_cmp_sd(__m128d, __m128d, const int);





extern int __cdecl _mm_comi_sd(__m128d, __m128d, const int);












extern __m128 __cdecl _mm_cmp_ss(__m128, __m128, const int);





extern int __cdecl _mm_comi_ss(__m128, __m128, const int);








extern __m256d __cdecl _mm256_cvtepi32_pd(__m128i);








extern __m256  __cdecl _mm256_cvtepi32_ps(__m256i);









extern __m128  __cdecl _mm256_cvtpd_ps(__m256d);








extern __m256i __cdecl _mm256_cvtps_epi32(__m256);









extern __m256d __cdecl _mm256_cvtps_pd(__m128);












extern __m128i __cdecl _mm256_cvttpd_epi32(__m256d);








extern __m128i __cdecl _mm256_cvtpd_epi32(__m256d);












extern __m256i __cdecl _mm256_cvttps_epi32(__m256);







extern __m128  __cdecl _mm256_extractf128_ps(__m256, const int);
extern __m128d __cdecl _mm256_extractf128_pd(__m256d, const int);
extern __m128i __cdecl _mm256_extractf128_si256(__m256i, const int);






extern void __cdecl _mm256_zeroall(void);







extern void __cdecl _mm256_zeroupper(void);









extern __m256  __cdecl _mm256_permutevar_ps(__m256, __m256i);
extern __m128  __cdecl _mm_permutevar_ps(__m128, __m128i);









extern __m256  __cdecl _mm256_permute_ps(__m256, int);
extern __m128  __cdecl _mm_permute_ps(__m128, int);









extern __m256d __cdecl _mm256_permutevar_pd(__m256d, __m256i);
extern __m128d __cdecl _mm_permutevar_pd(__m128d, __m128i);









extern __m256d __cdecl _mm256_permute_pd(__m256d, int);
extern __m128d __cdecl _mm_permute_pd(__m128d, int);








extern __m256  __cdecl _mm256_permute2f128_ps(__m256, __m256, int);
extern __m256d __cdecl _mm256_permute2f128_pd(__m256d, __m256d, int);
extern __m256i __cdecl _mm256_permute2f128_si256(__m256i, __m256i, int);








extern __m256  __cdecl _mm256_broadcast_ss(float const *);
extern __m128  __cdecl _mm_broadcast_ss(float const *);







extern __m256d __cdecl _mm256_broadcast_sd(double const *);







extern __m256  __cdecl _mm256_broadcast_ps(__m128 const *);
extern __m256d __cdecl _mm256_broadcast_pd(__m128d const *);









extern __m256  __cdecl _mm256_insertf128_ps(__m256, __m128, int);
extern __m256d __cdecl _mm256_insertf128_pd(__m256d, __m128d, int);
extern __m256i __cdecl _mm256_insertf128_si256(__m256i, __m128i, int);








extern __m256d __cdecl _mm256_load_pd(double const *);
extern void    __cdecl _mm256_store_pd(double *, __m256d);








extern __m256  __cdecl _mm256_load_ps(float const *);
extern void    __cdecl _mm256_store_ps(float *, __m256);








extern __m256d __cdecl _mm256_loadu_pd(double const *);
extern void    __cdecl _mm256_storeu_pd(double *, __m256d);








extern __m256  __cdecl _mm256_loadu_ps(float const *);
extern void    __cdecl _mm256_storeu_ps(float *, __m256);








extern __m256i __cdecl _mm256_load_si256(__m256i const *);
extern void    __cdecl _mm256_store_si256(__m256i *, __m256i);








extern __m256i __cdecl _mm256_loadu_si256(__m256i const *);
extern void    __cdecl _mm256_storeu_si256(__m256i *, __m256i);







































































extern __m256d __cdecl _mm256_maskload_pd(double const *, __m256i);
extern void    __cdecl _mm256_maskstore_pd(double *, __m256i, __m256d);
extern __m128d __cdecl _mm_maskload_pd(double const *, __m128i);
extern void    __cdecl _mm_maskstore_pd(double *, __m128i, __m128d);



















extern __m256  __cdecl _mm256_maskload_ps(float const *, __m256i);
extern void    __cdecl _mm256_maskstore_ps(float *, __m256i, __m256);
extern __m128  __cdecl _mm_maskload_ps(float const *, __m128i);
extern void    __cdecl _mm_maskstore_ps(float *, __m128i, __m128);







extern __m256  __cdecl _mm256_movehdup_ps(__m256);







extern __m256  __cdecl _mm256_moveldup_ps(__m256);







extern __m256d __cdecl _mm256_movedup_pd(__m256d);









extern __m256i __cdecl _mm256_lddqu_si256(__m256i const *);







extern void    __cdecl _mm256_stream_si256(__m256i *, __m256i);








extern void    __cdecl _mm256_stream_pd(double *, __m256d);








extern void    __cdecl _mm256_stream_ps(float *, __m256);









extern __m256  __cdecl _mm256_rcp_ps(__m256);










extern __m256  __cdecl _mm256_rsqrt_ps(__m256);








extern __m256d __cdecl _mm256_sqrt_pd(__m256d);








extern __m256  __cdecl _mm256_sqrt_ps(__m256);












extern __m256d __cdecl _mm256_round_pd(__m256d, int);














extern __m256  __cdecl _mm256_round_ps(__m256, int);









extern __m256d __cdecl _mm256_unpackhi_pd(__m256d, __m256d);







extern __m256  __cdecl _mm256_unpackhi_ps(__m256, __m256);







extern __m256d __cdecl _mm256_unpacklo_pd(__m256d, __m256d);







extern __m256  __cdecl _mm256_unpacklo_ps(__m256, __m256);









extern int     __cdecl _mm256_testz_si256(__m256i, __m256i);



extern int     __cdecl _mm256_testc_si256(__m256i, __m256i);



extern int     __cdecl _mm256_testnzc_si256(__m256i, __m256i);














extern int     __cdecl _mm256_testz_pd(__m256d, __m256d);
extern int     __cdecl _mm256_testc_pd(__m256d, __m256d);
extern int     __cdecl _mm256_testnzc_pd(__m256d, __m256d);
extern int     __cdecl _mm_testz_pd(__m128d, __m128d);
extern int     __cdecl _mm_testc_pd(__m128d, __m128d);
extern int     __cdecl _mm_testnzc_pd(__m128d, __m128d);












extern int     __cdecl _mm256_testz_ps(__m256, __m256);
extern int     __cdecl _mm256_testc_ps(__m256, __m256);
extern int     __cdecl _mm256_testnzc_ps(__m256, __m256);
extern int     __cdecl _mm_testz_ps(__m128, __m128);
extern int     __cdecl _mm_testc_ps(__m128, __m128);
extern int     __cdecl _mm_testnzc_ps(__m128, __m128);








extern int     __cdecl _mm256_movemask_pd(__m256d);








extern int     __cdecl _mm256_movemask_ps(__m256);




extern __m256d __cdecl _mm256_setzero_pd(void);
extern __m256  __cdecl _mm256_setzero_ps(void);
extern __m256i __cdecl _mm256_setzero_si256(void);




extern __m256d __cdecl _mm256_set_pd(double, double, double, double);
extern __m256  __cdecl _mm256_set_ps(float, float, float, float,
                                            float, float, float, float);
extern __m256i __cdecl _mm256_set_epi8(char, char, char, char,
                                              char, char, char, char,
                                              char, char, char, char,
                                              char, char, char, char,
                                              char, char, char, char,
                                              char, char, char, char,
                                              char, char, char, char,
                                              char, char, char, char);
extern __m256i __cdecl _mm256_set_epi16(short, short, short, short,
                                               short, short, short, short,
                                               short, short, short, short,
                                               short, short, short, short);
extern __m256i __cdecl _mm256_set_epi32(int, int, int, int,
                                               int, int, int, int);
extern __m256i __cdecl _mm256_set_epi64x(__int64, __int64,
                                                __int64, __int64);










extern __m256d __cdecl _mm256_setr_pd(double, double, double, double);
extern __m256  __cdecl _mm256_setr_ps(float, float, float, float,
                                             float, float, float, float);
extern __m256i __cdecl _mm256_setr_epi8(char, char, char, char,
                                               char, char, char, char,
                                               char, char, char, char,
                                               char, char, char, char,
                                               char, char, char, char,
                                               char, char, char, char,
                                               char, char, char, char,
                                               char, char, char, char);
extern __m256i __cdecl _mm256_setr_epi16(short, short, short, short,
                                                short, short, short, short,
                                                short, short, short, short,
                                                short, short, short, short);
extern __m256i __cdecl _mm256_setr_epi32(int, int, int, int,
                                                int, int, int, int);
extern __m256i __cdecl _mm256_setr_epi64x(__int64, __int64,
                                                 __int64, __int64);







extern __m256d __cdecl _mm256_set1_pd(double);
extern __m256  __cdecl _mm256_set1_ps(float);
extern __m256i __cdecl _mm256_set1_epi8(char);
extern __m256i __cdecl _mm256_set1_epi16(short);
extern __m256i __cdecl _mm256_set1_epi32(int);
extern __m256i __cdecl _mm256_set1_epi64x(long long);







extern __m256  __cdecl _mm256_castpd_ps(__m256d);
extern __m256d __cdecl _mm256_castps_pd(__m256);
extern __m256i __cdecl _mm256_castps_si256(__m256);
extern __m256i __cdecl _mm256_castpd_si256(__m256d);
extern __m256  __cdecl _mm256_castsi256_ps(__m256i);
extern __m256d __cdecl _mm256_castsi256_pd(__m256i);
extern __m128  __cdecl _mm256_castps256_ps128(__m256);
extern __m128d __cdecl _mm256_castpd256_pd128(__m256d);
extern __m128i __cdecl _mm256_castsi256_si128(__m256i);
extern __m256  __cdecl _mm256_castps128_ps256(__m128);
extern __m256d __cdecl _mm256_castpd128_pd256(__m128d);
extern __m256i __cdecl _mm256_castsi128_si256(__m128i);






extern __m128  __cdecl _mm_cvtph_ps(__m128i);
extern __m256  __cdecl _mm256_cvtph_ps(__m128i);
extern __m128i __cdecl _mm_cvtps_ph(__m128 , const int );
extern __m128i __cdecl _mm256_cvtps_ph(__m256, int);




















extern unsigned __int64 __cdecl _xgetbv(unsigned int);


extern void __cdecl _xsetbv(unsigned int, unsigned __int64);






extern void __cdecl _xsave(void *, unsigned __int64);

extern void __cdecl _xsave64(void *, unsigned __int64);







extern void __cdecl _xsaveopt(void *, unsigned __int64);

extern void __cdecl _xsaveopt64(void *, unsigned __int64);






extern void __cdecl _xsavec(void *, unsigned __int64);

extern void __cdecl _xsavec64(void *, unsigned __int64);







extern void __cdecl _xrstor(void const *, unsigned __int64);

extern void __cdecl _xrstor64(void const *, unsigned __int64);







extern void __cdecl _xsaves(void *, unsigned __int64);

extern void __cdecl _xsaves64(void *, unsigned __int64);







extern void __cdecl _xrstors(void const *, unsigned __int64);

extern void __cdecl _xrstors64(void const *, unsigned __int64);






extern void __cdecl _fxsave(void *);

extern void __cdecl _fxsave64(void *);






extern void __cdecl _fxrstor(void const *);

extern void __cdecl _fxrstor64(void const *);








extern int __cdecl _rdrand16_step(unsigned short *);
extern int __cdecl _rdrand32_step(unsigned int *);

extern int __cdecl _rdrand64_step(unsigned __int64 *);






extern unsigned int     __cdecl _readfsbase_u32();
extern unsigned int     __cdecl _readgsbase_u32();
extern unsigned __int64 __cdecl _readfsbase_u64();
extern unsigned __int64 __cdecl _readgsbase_u64();




extern void __cdecl _writefsbase_u32(unsigned int);
extern void __cdecl _writegsbase_u32(unsigned int);
extern void __cdecl _writefsbase_u64(unsigned __int64);
extern void __cdecl _writegsbase_u64(unsigned __int64);





extern __m128  __cdecl _mm_fmadd_ps(__m128, __m128, __m128);
extern __m128d __cdecl _mm_fmadd_pd(__m128d, __m128d, __m128d);
extern __m128  __cdecl _mm_fmadd_ss(__m128, __m128, __m128);
extern __m128d __cdecl _mm_fmadd_sd(__m128d, __m128d, __m128d);
extern __m128  __cdecl _mm_fmsub_ps(__m128, __m128, __m128);
extern __m128d __cdecl _mm_fmsub_pd(__m128d, __m128d, __m128d);
extern __m128  __cdecl _mm_fmsub_ss(__m128, __m128, __m128);
extern __m128d __cdecl _mm_fmsub_sd(__m128d, __m128d, __m128d);
extern __m128  __cdecl _mm_fnmadd_ps(__m128, __m128, __m128);
extern __m128d __cdecl _mm_fnmadd_pd(__m128d, __m128d, __m128d);
extern __m128  __cdecl _mm_fnmadd_ss(__m128, __m128, __m128);
extern __m128d __cdecl _mm_fnmadd_sd(__m128d, __m128d, __m128d);
extern __m128  __cdecl _mm_fnmsub_ps(__m128, __m128, __m128);
extern __m128d __cdecl _mm_fnmsub_pd(__m128d, __m128d, __m128d);
extern __m128  __cdecl _mm_fnmsub_ss(__m128, __m128, __m128);
extern __m128d __cdecl _mm_fnmsub_sd(__m128d, __m128d, __m128d);

extern __m256  __cdecl _mm256_fmadd_ps(__m256, __m256, __m256);
extern __m256d __cdecl _mm256_fmadd_pd(__m256d, __m256d, __m256d);
extern __m256  __cdecl _mm256_fmsub_ps(__m256, __m256, __m256);
extern __m256d __cdecl _mm256_fmsub_pd(__m256d, __m256d, __m256d);
extern __m256  __cdecl _mm256_fnmadd_ps(__m256, __m256, __m256);
extern __m256d __cdecl _mm256_fnmadd_pd(__m256d, __m256d, __m256d);
extern __m256  __cdecl _mm256_fnmsub_ps(__m256, __m256, __m256);
extern __m256d __cdecl _mm256_fnmsub_pd(__m256d, __m256d, __m256d);





extern __m128  __cdecl _mm_fmaddsub_ps(__m128, __m128, __m128);
extern __m128d __cdecl _mm_fmaddsub_pd(__m128d, __m128d, __m128d);
extern __m128  __cdecl _mm_fmsubadd_ps(__m128, __m128, __m128);
extern __m128d __cdecl _mm_fmsubadd_pd(__m128d, __m128d, __m128d);

extern __m256  __cdecl _mm256_fmaddsub_ps(__m256, __m256, __m256);
extern __m256d __cdecl _mm256_fmaddsub_pd(__m256d, __m256d, __m256d);
extern __m256  __cdecl _mm256_fmsubadd_ps(__m256, __m256, __m256);
extern __m256d __cdecl _mm256_fmsubadd_pd(__m256d, __m256d, __m256d);





extern __m256i __cdecl _mm256_cmpeq_epi8(__m256i, __m256i);
extern __m256i __cdecl _mm256_cmpeq_epi16(__m256i, __m256i);
extern __m256i __cdecl _mm256_cmpeq_epi32(__m256i, __m256i);
extern __m256i __cdecl _mm256_cmpeq_epi64(__m256i, __m256i);

extern __m256i __cdecl _mm256_cmpgt_epi8(__m256i, __m256i);
extern __m256i __cdecl _mm256_cmpgt_epi16(__m256i, __m256i);
extern __m256i __cdecl _mm256_cmpgt_epi32(__m256i, __m256i);
extern __m256i __cdecl _mm256_cmpgt_epi64(__m256i, __m256i);





extern __m256i __cdecl _mm256_max_epi8(__m256i, __m256i);
extern __m256i __cdecl _mm256_max_epi16(__m256i, __m256i);
extern __m256i __cdecl _mm256_max_epi32(__m256i, __m256i);
extern __m256i __cdecl _mm256_max_epu8(__m256i, __m256i);
extern __m256i __cdecl _mm256_max_epu16(__m256i, __m256i);
extern __m256i __cdecl _mm256_max_epu32(__m256i, __m256i);

extern __m256i __cdecl _mm256_min_epi8(__m256i, __m256i);
extern __m256i __cdecl _mm256_min_epi16(__m256i, __m256i);
extern __m256i __cdecl _mm256_min_epi32(__m256i, __m256i);
extern __m256i __cdecl _mm256_min_epu8(__m256i, __m256i);
extern __m256i __cdecl _mm256_min_epu16(__m256i, __m256i);
extern __m256i __cdecl _mm256_min_epu32(__m256i, __m256i);





extern __m256i __cdecl _mm256_and_si256(__m256i, __m256i);
extern __m256i __cdecl _mm256_andnot_si256(__m256i, __m256i);
extern __m256i __cdecl _mm256_or_si256(__m256i, __m256i);
extern __m256i __cdecl _mm256_xor_si256(__m256i, __m256i);





extern __m256i __cdecl _mm256_abs_epi8(__m256i);
extern __m256i __cdecl _mm256_abs_epi16(__m256i);
extern __m256i __cdecl _mm256_abs_epi32(__m256i);

extern __m256i __cdecl _mm256_add_epi8(__m256i, __m256i);
extern __m256i __cdecl _mm256_add_epi16(__m256i, __m256i);
extern __m256i __cdecl _mm256_add_epi32(__m256i, __m256i);
extern __m256i __cdecl _mm256_add_epi64(__m256i, __m256i);

extern __m256i __cdecl _mm256_adds_epi8(__m256i, __m256i);
extern __m256i __cdecl _mm256_adds_epi16(__m256i, __m256i);
extern __m256i __cdecl _mm256_adds_epu8(__m256i, __m256i);
extern __m256i __cdecl _mm256_adds_epu16(__m256i, __m256i);

extern __m256i __cdecl _mm256_sub_epi8(__m256i, __m256i);
extern __m256i __cdecl _mm256_sub_epi16(__m256i, __m256i);
extern __m256i __cdecl _mm256_sub_epi32(__m256i, __m256i);
extern __m256i __cdecl _mm256_sub_epi64(__m256i, __m256i);

extern __m256i __cdecl _mm256_subs_epi8(__m256i, __m256i);
extern __m256i __cdecl _mm256_subs_epi16(__m256i, __m256i);
extern __m256i __cdecl _mm256_subs_epu8(__m256i, __m256i);
extern __m256i __cdecl _mm256_subs_epu16(__m256i, __m256i);

extern __m256i __cdecl _mm256_avg_epu8(__m256i, __m256i);
extern __m256i __cdecl _mm256_avg_epu16(__m256i, __m256i);

extern __m256i __cdecl _mm256_hadd_epi16(__m256i, __m256i);
extern __m256i __cdecl _mm256_hadd_epi32(__m256i, __m256i);
extern __m256i __cdecl _mm256_hadds_epi16(__m256i, __m256i);

extern __m256i __cdecl _mm256_hsub_epi16(__m256i, __m256i);
extern __m256i __cdecl _mm256_hsub_epi32(__m256i, __m256i);
extern __m256i __cdecl _mm256_hsubs_epi16(__m256i, __m256i);

extern __m256i __cdecl _mm256_madd_epi16(__m256i, __m256i);
extern __m256i __cdecl _mm256_maddubs_epi16(__m256i, __m256i);

extern __m256i __cdecl _mm256_mulhi_epi16(__m256i, __m256i);
extern __m256i __cdecl _mm256_mulhi_epu16(__m256i, __m256i);

extern __m256i __cdecl _mm256_mullo_epi16(__m256i, __m256i);
extern __m256i __cdecl _mm256_mullo_epi32(__m256i, __m256i);

extern __m256i __cdecl _mm256_mul_epu32(__m256i, __m256i);
extern __m256i __cdecl _mm256_mul_epi32(__m256i, __m256i);

extern __m256i __cdecl _mm256_sign_epi8(__m256i, __m256i);
extern __m256i __cdecl _mm256_sign_epi16(__m256i, __m256i);
extern __m256i __cdecl _mm256_sign_epi32(__m256i, __m256i);

extern __m256i __cdecl _mm256_mulhrs_epi16(__m256i, __m256i);

extern __m256i __cdecl _mm256_sad_epu8(__m256i, __m256i);
extern __m256i __cdecl _mm256_mpsadbw_epu8(__m256i, __m256i, const int);





extern __m256i __cdecl _mm256_slli_si256(__m256i, const int);

extern __m256i __cdecl _mm256_srli_si256(__m256i, const int);


extern __m256i __cdecl _mm256_sll_epi16(__m256i, __m128i);
extern __m256i __cdecl _mm256_sll_epi32(__m256i, __m128i);
extern __m256i __cdecl _mm256_sll_epi64(__m256i, __m128i);

extern __m256i __cdecl _mm256_slli_epi16(__m256i, int);
extern __m256i __cdecl _mm256_slli_epi32(__m256i, int);
extern __m256i __cdecl _mm256_slli_epi64(__m256i, int);

extern __m256i __cdecl _mm256_sllv_epi32(__m256i, __m256i);
extern __m256i __cdecl _mm256_sllv_epi64(__m256i, __m256i);

extern __m128i __cdecl _mm_sllv_epi32(__m128i, __m128i);
extern __m128i __cdecl _mm_sllv_epi64(__m128i, __m128i);

extern __m256i __cdecl _mm256_sra_epi16(__m256i, __m128i);
extern __m256i __cdecl _mm256_sra_epi32(__m256i, __m128i);

extern __m256i __cdecl _mm256_srai_epi16(__m256i, int);
extern __m256i __cdecl _mm256_srai_epi32(__m256i, int);

extern __m256i __cdecl _mm256_srav_epi32(__m256i, __m256i);

extern __m128i __cdecl _mm_srav_epi32(__m128i, __m128i);

extern __m256i __cdecl _mm256_srl_epi16(__m256i, __m128i);
extern __m256i __cdecl _mm256_srl_epi32(__m256i, __m128i);
extern __m256i __cdecl _mm256_srl_epi64(__m256i, __m128i);

extern __m256i __cdecl _mm256_srli_epi16(__m256i, int);
extern __m256i __cdecl _mm256_srli_epi32(__m256i, int);
extern __m256i __cdecl _mm256_srli_epi64(__m256i, int);

extern __m256i __cdecl _mm256_srlv_epi32(__m256i, __m256i);
extern __m256i __cdecl _mm256_srlv_epi64(__m256i, __m256i);

extern __m128i __cdecl _mm_srlv_epi32(__m128i, __m128i);
extern __m128i __cdecl _mm_srlv_epi64(__m128i, __m128i);





extern __m128i __cdecl _mm_blend_epi32(__m128i, __m128i, const int);

extern __m256i __cdecl _mm256_blend_epi32(__m256i,__m256i, const int);

extern __m256i __cdecl _mm256_alignr_epi8(__m256i, __m256i, const int);

extern __m256i __cdecl _mm256_blendv_epi8(__m256i, __m256i, __m256i);
extern __m256i __cdecl _mm256_blend_epi16(__m256i, __m256i, const int);

extern __m256i __cdecl _mm256_packs_epi16(__m256i, __m256i);
extern __m256i __cdecl _mm256_packs_epi32(__m256i, __m256i);
extern __m256i __cdecl _mm256_packus_epi16(__m256i, __m256i);
extern __m256i __cdecl _mm256_packus_epi32(__m256i, __m256i);

extern __m256i __cdecl _mm256_unpackhi_epi8(__m256i, __m256i);
extern __m256i __cdecl _mm256_unpackhi_epi16(__m256i, __m256i);
extern __m256i __cdecl _mm256_unpackhi_epi32(__m256i, __m256i);
extern __m256i __cdecl _mm256_unpackhi_epi64(__m256i, __m256i);

extern __m256i __cdecl _mm256_unpacklo_epi8(__m256i, __m256i);
extern __m256i __cdecl _mm256_unpacklo_epi16(__m256i, __m256i);
extern __m256i __cdecl _mm256_unpacklo_epi32(__m256i, __m256i);
extern __m256i __cdecl _mm256_unpacklo_epi64(__m256i, __m256i);

extern __m256i __cdecl _mm256_shuffle_epi8(__m256i, __m256i);
extern __m256i __cdecl _mm256_shuffle_epi32(__m256i, const int);

extern __m256i __cdecl _mm256_shufflehi_epi16(__m256i, const int);
extern __m256i __cdecl _mm256_shufflelo_epi16(__m256i, const int);

extern __m128i __cdecl _mm256_extracti128_si256(__m256i, const int);
extern __m256i __cdecl _mm256_inserti128_si256(__m256i, __m128i, const int);





extern __m128  __cdecl _mm_broadcastss_ps(__m128);
extern __m128d __cdecl _mm_broadcastsd_pd(__m128d);

extern __m128i __cdecl _mm_broadcastb_epi8(__m128i);
extern __m128i __cdecl _mm_broadcastw_epi16(__m128i);
extern __m128i __cdecl _mm_broadcastd_epi32(__m128i);
extern __m128i __cdecl _mm_broadcastq_epi64(__m128i);

extern __m256  __cdecl _mm256_broadcastss_ps(__m128);
extern __m256d __cdecl _mm256_broadcastsd_pd(__m128d);

extern __m256i __cdecl _mm256_broadcastb_epi8(__m128i);
extern __m256i __cdecl _mm256_broadcastw_epi16(__m128i);
extern __m256i __cdecl _mm256_broadcastd_epi32(__m128i);
extern __m256i __cdecl _mm256_broadcastq_epi64(__m128i);

extern __m256i __cdecl _mm256_broadcastsi128_si256(__m128i);






extern __m256i __cdecl _mm256_cvtepi8_epi16(__m128i);
extern __m256i __cdecl _mm256_cvtepi8_epi32(__m128i);
extern __m256i __cdecl _mm256_cvtepi8_epi64(__m128i);
extern __m256i __cdecl _mm256_cvtepi16_epi32(__m128i);
extern __m256i __cdecl _mm256_cvtepi16_epi64(__m128i);
extern __m256i __cdecl _mm256_cvtepi32_epi64(__m128i);

extern __m256i __cdecl _mm256_cvtepu8_epi16(__m128i);
extern __m256i __cdecl _mm256_cvtepu8_epi32(__m128i);
extern __m256i __cdecl _mm256_cvtepu8_epi64(__m128i);
extern __m256i __cdecl _mm256_cvtepu16_epi32(__m128i);
extern __m256i __cdecl _mm256_cvtepu16_epi64(__m128i);
extern __m256i __cdecl _mm256_cvtepu32_epi64(__m128i);






extern int __cdecl _mm256_movemask_epi8(__m256i);





extern __m128i __cdecl _mm_maskload_epi32(int const * ,
                                          __m128i     );
extern __m128i __cdecl _mm_maskload_epi64(__int64 const * ,
                                          __m128i         );

extern void __cdecl _mm_maskstore_epi32(int *   ,
                                        __m128i ,
                                        __m128i );
extern void __cdecl _mm_maskstore_epi64(__int64 * ,
                                        __m128i   ,
                                        __m128i   );

extern __m256i __cdecl _mm256_maskload_epi32(int const * ,
                                             __m256i     );
extern __m256i __cdecl _mm256_maskload_epi64(__int64 const * ,
                                             __m256i         );

extern void __cdecl _mm256_maskstore_epi32(int *   ,
                                           __m256i ,
                                           __m256i );
extern void __cdecl _mm256_maskstore_epi64(__int64 * ,
                                           __m256i   ,
                                           __m256i   );





extern __m256i __cdecl _mm256_permutevar8x32_epi32(__m256i, __m256i);
extern __m256  __cdecl _mm256_permutevar8x32_ps(__m256, __m256i);

extern __m256i __cdecl _mm256_permute4x64_epi64(__m256i, const int);
extern __m256d __cdecl _mm256_permute4x64_pd(__m256d, const int);

extern __m256i __cdecl _mm256_permute2x128_si256(__m256i, __m256i, const int);





extern __m256i  __cdecl _mm256_stream_load_si256(__m256i const *);






extern __m256d __cdecl _mm256_mask_i32gather_pd(__m256d        ,
                                                double const * ,
                                                __m128i        ,
                                                __m256d        ,
                                                const int      );
extern __m256  __cdecl _mm256_mask_i32gather_ps(__m256         ,
                                                float const *  ,
                                                __m256i        ,
                                                __m256         ,
                                                const int      );
extern __m256d __cdecl _mm256_mask_i64gather_pd(__m256d        ,
                                                double const * ,
                                                __m256i        ,
                                                __m256d        ,
                                                const int      );
extern __m128  __cdecl _mm256_mask_i64gather_ps(__m128         ,
                                                float const *  ,
                                                __m256i        ,
                                                __m128         ,
                                                const int      );

extern __m128d __cdecl _mm_mask_i32gather_pd(__m128d        ,
                                             double const * ,
                                             __m128i        ,
                                             __m128d        ,
                                             const int      );
extern __m128  __cdecl _mm_mask_i32gather_ps(__m128         ,
                                             float const *  ,
                                             __m128i        ,
                                             __m128         ,
                                             const int      );
extern __m128d __cdecl _mm_mask_i64gather_pd(__m128d        ,
                                             double const * ,
                                             __m128i        ,
                                             __m128d        ,
                                             const int      );
extern __m128  __cdecl _mm_mask_i64gather_ps(__m128         ,
                                             float const *  ,
                                             __m128i        ,
                                             __m128         ,
                                             const int      );


extern __m256i __cdecl _mm256_mask_i32gather_epi32(__m256i     ,
                                                   int const * ,
                                                   __m256i     ,
                                                   __m256i     ,
                                                   const int   );
extern __m256i __cdecl _mm256_mask_i32gather_epi64(__m256i     ,
                                                   __int64 const * ,
                                                   __m128i     ,
                                                   __m256i     ,
                                                   const int   );
extern __m128i __cdecl _mm256_mask_i64gather_epi32(__m128i     ,
                                                   int     const * ,
                                                   __m256i     ,
                                                   __m128i     ,
                                                   const int   );
extern __m256i __cdecl _mm256_mask_i64gather_epi64(__m256i     ,
                                                   __int64 const * ,
                                                   __m256i     ,
                                                   __m256i     ,
                                                   const int   );

extern __m128i __cdecl _mm_mask_i32gather_epi32(__m128i         ,
                                                int const *     ,
                                                __m128i         ,
                                                __m128i         ,
                                                const int       );
extern __m128i __cdecl _mm_mask_i32gather_epi64(__m128i         ,
                                                __int64 const * ,
                                                __m128i         ,
                                                __m128i         ,
                                                const int       );
extern __m128i __cdecl _mm_mask_i64gather_epi32(__m128i         ,
                                                int     const * ,
                                                __m128i         ,
                                                __m128i         ,
                                                const int       );
extern __m128i __cdecl _mm_mask_i64gather_epi64(__m128i         ,
                                                __int64 const * ,
                                                __m128i         ,
                                                __m128i         ,
                                                const int       );





extern __m256d __cdecl _mm256_i32gather_pd(double const * ,
                                           __m128i        ,
                                           const int      );
extern __m256  __cdecl _mm256_i32gather_ps(float  const * ,
                                           __m256i        ,
                                           const int      );
extern __m256d __cdecl _mm256_i64gather_pd(double const * ,
                                           __m256i        ,
                                           const int      );
extern __m128  __cdecl _mm256_i64gather_ps(float  const * ,
                                           __m256i        ,
                                           const int      );

extern __m128d __cdecl _mm_i32gather_pd(double const * ,
                                        __m128i        ,
                                        const int      );
extern __m128  __cdecl _mm_i32gather_ps(float  const * ,
                                        __m128i        ,
                                        const int      );
extern __m128d __cdecl _mm_i64gather_pd(double const * ,
                                        __m128i        ,
                                        const int      );
extern __m128  __cdecl _mm_i64gather_ps(float  const * ,
                                        __m128i        ,
                                        const int      );

extern __m256i __cdecl _mm256_i32gather_epi32(int const *     ,
                                              __m256i         ,
                                              const int       );
extern __m256i __cdecl _mm256_i32gather_epi64(__int64 const * ,
                                              __m128i         ,
                                              const int       );
extern __m128i __cdecl _mm256_i64gather_epi32(int const *     ,
                                              __m256i         ,
                                              const int       );
extern __m256i __cdecl _mm256_i64gather_epi64(__int64 const * ,
                                              __m256i         ,
                                              const int       );

extern __m128i __cdecl _mm_i32gather_epi32(int const *     ,
                                           __m128i         ,
                                           const int       );
extern __m128i __cdecl _mm_i32gather_epi64(__int64 const * ,
                                           __m128i         ,
                                           const int       );
extern __m128i __cdecl _mm_i64gather_epi32(int     const * ,
                                           __m128i         ,
                                           const int       );
extern __m128i __cdecl _mm_i64gather_epi64(__int64 const * ,
                                           __m128i         ,
                                           const int       );







extern unsigned int     _bextr_u32(unsigned int ,
                                   unsigned int ,
                                   unsigned int );
extern unsigned int     _blsi_u32(unsigned int);
extern unsigned int     _blsmsk_u32(unsigned int);
extern unsigned int     _blsr_u32(unsigned int);
extern unsigned int     _bzhi_u32(unsigned int ,
                                  unsigned int );
extern unsigned int     _mulx_u32(unsigned int ,
                                  unsigned int ,
                                  unsigned int * );
extern unsigned int     _pdep_u32(unsigned int ,
                                  unsigned int );
extern unsigned int     _pext_u32(unsigned int ,
                                  unsigned int );
extern unsigned int     _rorx_u32(unsigned int ,
                                  const unsigned int );
extern int              _sarx_i32(int ,
                                  unsigned int );
extern unsigned int     _shlx_u32(unsigned int ,
                                  unsigned int );
extern unsigned int     _shrx_u32(unsigned int ,
                                          unsigned int );


extern unsigned __int64 _bextr_u64(unsigned __int64 ,
                                   unsigned int ,
                                   unsigned int );
extern unsigned __int64 _blsi_u64(unsigned __int64);
extern unsigned __int64 _blsmsk_u64(unsigned __int64);
extern unsigned __int64 _blsr_u64(unsigned __int64);
extern unsigned __int64 _bzhi_u64(unsigned __int64 ,
                                  unsigned int );
extern unsigned __int64 _mulx_u64(unsigned __int64 ,
                                  unsigned __int64 ,
                                  unsigned __int64 * );
extern unsigned __int64 _pdep_u64(unsigned __int64 ,
                                  unsigned __int64 );
extern unsigned __int64 _pext_u64(unsigned __int64 ,
                                  unsigned __int64 );
extern unsigned __int64 _rorx_u64(unsigned __int64 ,
                                  const unsigned int );
extern __int64          _sarx_i64(__int64 ,
                                  unsigned int );
extern unsigned __int64 _shlx_u64(unsigned __int64 ,
                                  unsigned int );
extern unsigned __int64 _shrx_u64(unsigned __int64 ,
                                          unsigned int );









extern unsigned int     _lzcnt_u32(unsigned int);

extern unsigned __int64 _lzcnt_u64(unsigned __int64);









extern unsigned int     _tzcnt_u32(unsigned int);

extern unsigned __int64 _tzcnt_u64(unsigned __int64);







extern void __cdecl _invpcid(unsigned int , void * );


extern void _Store_HLERelease(long volatile *,long);
extern void _StorePointer_HLERelease(void * volatile *,void *);

extern long _InterlockedExchange_HLEAcquire(long volatile *,long);
extern long _InterlockedExchange_HLERelease(long volatile *,long);
extern void * _InterlockedExchangePointer_HLEAcquire(void *volatile *,void *);
extern void * _InterlockedExchangePointer_HLERelease(void *volatile *,void *);

extern long _InterlockedCompareExchange_HLEAcquire(long volatile *,long,long);
extern long _InterlockedCompareExchange_HLERelease(long volatile *,long,long);
extern __int64 _InterlockedCompareExchange64_HLEAcquire(__int64 volatile *,__int64,__int64);
extern __int64 _InterlockedCompareExchange64_HLERelease(__int64 volatile *,__int64,__int64);
extern void * _InterlockedCompareExchangePointer_HLEAcquire(void *volatile *,void *,void *);
extern void * _InterlockedCompareExchangePointer_HLERelease(void *volatile *,void *,void *);

extern long _InterlockedExchangeAdd_HLEAcquire(long volatile *,long);
extern long _InterlockedExchangeAdd_HLERelease(long volatile *,long);

extern long _InterlockedAnd_HLEAcquire(long volatile *,long);
extern long _InterlockedAnd_HLERelease(long volatile *,long);
extern long _InterlockedOr_HLEAcquire(long volatile *,long);
extern long _InterlockedOr_HLERelease(long volatile *,long);
extern long _InterlockedXor_HLEAcquire(long volatile *,long);
extern long _InterlockedXor_HLERelease(long volatile *,long);

extern unsigned char _interlockedbittestandset_HLEAcquire(long *,long);
extern unsigned char _interlockedbittestandset_HLERelease(long *,long);
extern unsigned char _interlockedbittestandreset_HLEAcquire(long *,long);
extern unsigned char _interlockedbittestandreset_HLERelease(long *,long);


extern void _Store64_HLERelease(__int64 volatile *,__int64);
extern __int64 _InterlockedExchange64_HLEAcquire(__int64 volatile *,__int64);
extern __int64 _InterlockedExchange64_HLERelease(__int64 volatile *,__int64);

extern __int64 _InterlockedExchangeAdd64_HLEAcquire(__int64 volatile *,__int64);
extern __int64 _InterlockedExchangeAdd64_HLERelease(__int64 volatile *,__int64);

extern __int64 _InterlockedAnd64_HLEAcquire(__int64 volatile *,__int64);
extern __int64 _InterlockedAnd64_HLERelease(__int64 volatile *,__int64);
extern __int64 _InterlockedOr64_HLEAcquire(__int64 volatile *,__int64);
extern __int64 _InterlockedOr64_HLERelease(__int64 volatile *,__int64);
extern __int64 _InterlockedXor64_HLEAcquire(__int64 volatile *,__int64);
extern __int64 _InterlockedXor64_HLERelease(__int64 volatile *,__int64);

extern unsigned char _interlockedbittestandset64_HLEAcquire(__int64 *,__int64);
extern unsigned char _interlockedbittestandset64_HLERelease(__int64 *,__int64);
extern unsigned char _interlockedbittestandreset64_HLEAcquire(__int64 *,__int64);
extern unsigned char _interlockedbittestandreset64_HLERelease(__int64 *,__int64);












extern unsigned int     __cdecl _xbegin(void);
extern void             __cdecl _xend(void);
extern void             __cdecl _xabort(const unsigned int);
extern unsigned char    __cdecl _xtest(void);








extern int __cdecl _rdseed16_step(unsigned short *);
extern int __cdecl _rdseed32_step(unsigned int *);

extern int __cdecl _rdseed64_step(unsigned __int64 *);











extern unsigned char __cdecl _addcarryx_u32(unsigned char ,
                                                   unsigned int ,
                                                   unsigned int ,
                                                   unsigned int * );



extern unsigned char __cdecl _addcarryx_u64(unsigned char ,
                                                   unsigned __int64 ,
                                                   unsigned __int64 ,
                                                   unsigned __int64 * );






extern unsigned short   __cdecl _load_be_u16(void const*);
extern unsigned int     __cdecl _load_be_u32(void const*);
extern unsigned __int64 __cdecl _load_be_u64(void const*);







extern void __cdecl _store_be_u16(void *, unsigned short);
extern void __cdecl _store_be_u32(void *, unsigned int);
extern void __cdecl _store_be_u64(void *, unsigned __int64);







extern __m128i __cdecl _mm_sha1msg1_epu32(__m128i, __m128i);
extern __m128i __cdecl _mm_sha1msg2_epu32(__m128i, __m128i);
extern __m128i __cdecl _mm_sha1nexte_epu32(__m128i, __m128i);
extern __m128i __cdecl _mm_sha1rnds4_epu32(__m128i, __m128i, const int);

extern __m128i __cdecl _mm_sha256msg1_epu32(__m128i, __m128i);
extern __m128i __cdecl _mm_sha256msg2_epu32(__m128i, __m128i);
extern __m128i __cdecl _mm_sha256rnds2_epu32(__m128i, __m128i, __m128i);




extern void * __cdecl _bnd_set_ptr_bounds(const void *, size_t);
extern void * __cdecl _bnd_narrow_ptr_bounds(const void *, const void *, size_t);
extern void * __cdecl _bnd_copy_ptr_bounds(const void *, const void *);
extern void * __cdecl _bnd_init_ptr_bounds(const void *);
extern void __cdecl _bnd_store_ptr_bounds(const void **, const void *);
extern void __cdecl _bnd_chk_ptr_lbounds(const void *);
extern void __cdecl _bnd_chk_ptr_ubounds(const void *);
extern void __cdecl _bnd_chk_ptr_bounds(const void *, size_t);
extern void * __cdecl _bnd_load_ptr_bounds(const void **, const void *);
extern const void * __cdecl _bnd_get_ptr_lbound(const void *);
extern const void * __cdecl _bnd_get_ptr_ubound(const void *);


}; 







        












#pragma once


















extern "C" { 






























































































__m128 _mm_macc_ps(__m128, __m128, __m128);
__m128d _mm_macc_pd(__m128d, __m128d, __m128d);
__m128 _mm_macc_ss(__m128, __m128, __m128);
__m128d _mm_macc_sd(__m128d, __m128d, __m128d);
__m128 _mm_maddsub_ps(__m128, __m128, __m128);
__m128d _mm_maddsub_pd(__m128d, __m128d, __m128d);
__m128 _mm_msubadd_ps(__m128, __m128, __m128);
__m128d _mm_msubadd_pd(__m128d, __m128d, __m128d);
__m128 _mm_msub_ps(__m128, __m128, __m128);
__m128d _mm_msub_pd(__m128d, __m128d, __m128d);
__m128 _mm_msub_ss(__m128, __m128, __m128);
__m128d _mm_msub_sd(__m128d, __m128d, __m128d);
__m128 _mm_nmacc_ps(__m128, __m128, __m128);
__m128d _mm_nmacc_pd(__m128d, __m128d, __m128d);
__m128 _mm_nmacc_ss(__m128, __m128, __m128);
__m128d _mm_nmacc_sd(__m128d, __m128d, __m128d);
__m128 _mm_nmsub_ps(__m128, __m128, __m128);
__m128d _mm_nmsub_pd(__m128d, __m128d, __m128d);
__m128 _mm_nmsub_ss(__m128, __m128, __m128);
__m128d _mm_nmsub_sd(__m128d, __m128d, __m128d);


__m128i _mm_maccs_epi16(__m128i, __m128i, __m128i);
__m128i _mm_macc_epi16(__m128i, __m128i, __m128i);
__m128i _mm_maccsd_epi16(__m128i, __m128i, __m128i);
__m128i _mm_maccd_epi16(__m128i, __m128i, __m128i);
__m128i _mm_maccs_epi32(__m128i, __m128i, __m128i);
__m128i _mm_macc_epi32(__m128i, __m128i, __m128i);
__m128i _mm_maccslo_epi32(__m128i, __m128i, __m128i);
__m128i _mm_macclo_epi32(__m128i, __m128i, __m128i);
__m128i _mm_maccshi_epi32(__m128i, __m128i, __m128i);
__m128i _mm_macchi_epi32(__m128i, __m128i, __m128i);
__m128i _mm_maddsd_epi16(__m128i, __m128i, __m128i);
__m128i _mm_maddd_epi16(__m128i, __m128i, __m128i);


__m128i _mm_haddw_epi8(__m128i);
__m128i _mm_haddd_epi8(__m128i);
__m128i _mm_haddq_epi8(__m128i);
__m128i _mm_haddd_epi16(__m128i);
__m128i _mm_haddq_epi16(__m128i);
__m128i _mm_haddq_epi32(__m128i);
__m128i _mm_haddw_epu8(__m128i);
__m128i _mm_haddd_epu8(__m128i);
__m128i _mm_haddq_epu8(__m128i);
__m128i _mm_haddd_epu16(__m128i);
__m128i _mm_haddq_epu16(__m128i);
__m128i _mm_haddq_epu32(__m128i);
__m128i _mm_hsubw_epi8(__m128i);
__m128i _mm_hsubd_epi16(__m128i);
__m128i _mm_hsubq_epi32(__m128i);


__m128i _mm_cmov_si128(__m128i, __m128i, __m128i);
__m128i _mm_perm_epi8(__m128i, __m128i, __m128i);


__m128i _mm_rot_epi8(__m128i, __m128i);
__m128i _mm_rot_epi16(__m128i, __m128i);
__m128i _mm_rot_epi32(__m128i, __m128i);
__m128i _mm_rot_epi64(__m128i, __m128i);
__m128i _mm_roti_epi8(__m128i, int);
__m128i _mm_roti_epi16(__m128i, int);
__m128i _mm_roti_epi32(__m128i, int);
__m128i _mm_roti_epi64(__m128i, int);
__m128i _mm_shl_epi8(__m128i, __m128i);
__m128i _mm_shl_epi16(__m128i, __m128i);
__m128i _mm_shl_epi32(__m128i, __m128i);
__m128i _mm_shl_epi64(__m128i, __m128i);
__m128i _mm_sha_epi8(__m128i, __m128i);
__m128i _mm_sha_epi16(__m128i, __m128i);
__m128i _mm_sha_epi32(__m128i, __m128i);
__m128i _mm_sha_epi64(__m128i, __m128i);



__m128i _mm_com_epu8(__m128i, __m128i, int);
__m128i _mm_com_epu16(__m128i, __m128i, int);
__m128i _mm_com_epu32(__m128i, __m128i, int);
__m128i _mm_com_epu64(__m128i, __m128i, int);
__m128i _mm_com_epi8(__m128i, __m128i, int);
__m128i _mm_com_epi16(__m128i, __m128i, int);
__m128i _mm_com_epi32(__m128i, __m128i, int);
__m128i _mm_com_epi64(__m128i, __m128i, int);



__m128 _mm_frcz_ps(__m128);
__m128d _mm_frcz_pd(__m128d);
__m128 _mm_frcz_ss(__m128, __m128);
__m128d _mm_frcz_sd(__m128d, __m128d);








__m128 _mm_permute2_ps(__m128, __m128, __m128i, int);
__m128d _mm_permute2_pd(__m128d, __m128d, __m128i, int);



__m256 _mm256_macc_ps(__m256, __m256, __m256);
__m256d _mm256_macc_pd(__m256d, __m256d, __m256d);
__m256 _mm256_maddsub_ps(__m256, __m256, __m256);
__m256d _mm256_maddsub_pd(__m256d, __m256d, __m256d);
__m256 _mm256_msubadd_ps(__m256, __m256, __m256);
__m256d _mm256_msubadd_pd(__m256d, __m256d, __m256d);
__m256 _mm256_msub_ps(__m256, __m256, __m256);
__m256d _mm256_msub_pd(__m256d, __m256d, __m256d);
__m256 _mm256_nmacc_ps(__m256, __m256, __m256);
__m256d _mm256_nmacc_pd(__m256d, __m256d, __m256d);
__m256 _mm256_nmsub_ps(__m256, __m256, __m256);
__m256d _mm256_nmsub_pd(__m256d, __m256d, __m256d);
__m256i _mm256_cmov_si256(__m256i, __m256i, __m256i);
__m256 _mm256_frcz_ps(__m256);
__m256d _mm256_frcz_pd(__m256d);
__m256 _mm256_permute2_ps(__m256, __m256, __m256i, int);
__m256d _mm256_permute2_pd(__m256d, __m256d, __m256i, int);


void __llwpcb(void *);
void *__slwpcb();
void __lwpval32(unsigned int, unsigned int, unsigned int);
unsigned char __lwpins32(unsigned int, unsigned int, unsigned int);

void __lwpval64(unsigned __int64, unsigned int, unsigned int);
unsigned char __lwpins64(unsigned __int64, unsigned int, unsigned int);



unsigned int _bextr_u32(unsigned int, unsigned int, unsigned int);
unsigned int _andn_u32(unsigned int, unsigned int);
unsigned int _tzcnt_u32(unsigned int);
unsigned int _lzcnt_u32(unsigned int);
unsigned int _blsr_u32(unsigned int);
unsigned int _blsmsk_u32(unsigned int);
unsigned int _blsi_u32(unsigned int);

unsigned __int64 _bextr_u64(unsigned __int64, unsigned int, unsigned int);
unsigned __int64 _andn_u64(unsigned __int64, unsigned __int64);
unsigned __int64 _tzcnt_u64(unsigned __int64);
unsigned __int64 _lzcnt_u64(unsigned __int64);
unsigned __int64 _blsr_u64(unsigned __int64);
unsigned __int64 _blsmsk_u64(unsigned __int64);
unsigned __int64 _blsi_u64(unsigned __int64);



unsigned int _bextri_u32(unsigned int, unsigned int);
unsigned int _blcfill_u32(unsigned int);
unsigned int _blsfill_u32(unsigned int);
unsigned int _blcs_u32(unsigned int);
unsigned int _tzmsk_u32(unsigned int);
unsigned int _blcic_u32(unsigned int);
unsigned int _blsic_u32(unsigned int);
unsigned int _t1mskc_u32(unsigned int);
unsigned int _blcmsk_u32(unsigned int);
unsigned int _blci_u32(unsigned int);

unsigned __int64 _bextri_u64(unsigned __int64, unsigned int);
unsigned __int64 _blcfill_u64(unsigned __int64);
unsigned __int64 _blsfill_u64(unsigned __int64);
unsigned __int64 _blcs_u64(unsigned __int64);
unsigned __int64 _tzmsk_u64(unsigned __int64);
unsigned __int64 _blcic_u64(unsigned __int64);
unsigned __int64 _blsic_u64(unsigned __int64);
unsigned __int64 _t1mskc_u64(unsigned __int64);
unsigned __int64 _blcmsk_u64(unsigned __int64);
unsigned __int64 _blci_u64(unsigned __int64);


void _mm_monitorx(void const *, unsigned int, unsigned int);
void _mm_mwaitx(unsigned int, unsigned int, unsigned int);

void _mm_clzero(void const *);


}; 






    

    



    




    






extern "C" {




















































































void * _AddressOfReturnAddress(void);
unsigned char _BitScanForward(unsigned long * _Index, unsigned long _Mask);
unsigned char _BitScanForward64(unsigned long * _Index, unsigned __int64 _Mask);

unsigned char _BitScanReverse(unsigned long * _Index, unsigned long _Mask);
unsigned char _BitScanReverse64(unsigned long * _Index, unsigned __int64 _Mask);
























long _InterlockedAnd(long volatile * _Value, long _Mask);
short _InterlockedAnd16(short volatile * _Value, short _Mask);


short _InterlockedAnd16_np(short volatile * _Value, short _Mask);

__int64 _InterlockedAnd64(__int64 volatile * _Value, __int64 _Mask);


__int64 _InterlockedAnd64_np(__int64 volatile * _Value, __int64 _Mask);

char _InterlockedAnd8(char volatile * _Value, char _Mask);


char _InterlockedAnd8_np(char volatile * _Value, char _Mask);



long _InterlockedAnd_np(long volatile * _Value, long _Mask);

long  _InterlockedCompareExchange(long volatile * _Destination, long _Exchange, long _Comparand);

unsigned char _InterlockedCompareExchange128(__int64 volatile * _Destination, __int64 _ExchangeHigh, __int64 _ExchangeLow, __int64 * _ComparandResult);


unsigned char _InterlockedCompareExchange128_np(__int64 volatile * _Destination, __int64 _ExchangeHigh, __int64 _ExchangeLow, __int64 * _ComparandResult);

short _InterlockedCompareExchange16(short volatile * _Destination, short _Exchange, short _Comparand);


short _InterlockedCompareExchange16_np(short volatile * _Destination, short _Exchange, short _Comparand);

__int64 _InterlockedCompareExchange64(__int64 volatile * _Destination, __int64 _Exchange, __int64 _Comparand);


__int64 _InterlockedCompareExchange64_np(__int64 volatile * _Destination, __int64 _Exchange, __int64 _Comparand);

char _InterlockedCompareExchange8(char volatile * _Destination, char _Exchange, char _Comparand);



void * _InterlockedCompareExchangePointer(void * volatile * _Destination, void * _Exchange, void * _Comparand);


void * _InterlockedCompareExchangePointer_np(void * volatile * _Destination, void * _Exchange, void * _Comparand);



long _InterlockedCompareExchange_np(long volatile * _Destination, long _Exchange, long _Comparand);

long  _InterlockedDecrement(long volatile * _Addend);

short _InterlockedDecrement16(short volatile * _Addend);



__int64 _InterlockedDecrement64(__int64 volatile * _Addend);






long  _InterlockedExchange(long volatile * _Target, long _Value);

short _InterlockedExchange16(short volatile * _Target, short _Value);



__int64 _InterlockedExchange64(__int64 volatile * _Target, __int64 _Value);



char _InterlockedExchange8(char volatile * _Target, char _Value);



long  _InterlockedExchangeAdd(long volatile * _Addend, long _Value);
short _InterlockedExchangeAdd16(short volatile * _Addend, short _Value);



__int64 _InterlockedExchangeAdd64(__int64 volatile * _Addend, __int64 _Value);



char _InterlockedExchangeAdd8(char volatile * _Addend, char _Value);






void * _InterlockedExchangePointer(void * volatile * _Target, void * _Value);






long  _InterlockedIncrement(long volatile * _Addend);

short _InterlockedIncrement16(short volatile * _Addend);



__int64 _InterlockedIncrement64(__int64 volatile * _Addend);






long _InterlockedOr(long volatile * _Value, long _Mask);
short _InterlockedOr16(short volatile * _Value, short _Mask);


short _InterlockedOr16_np(short volatile * _Value, short _Mask);

__int64 _InterlockedOr64(__int64 volatile * _Value, __int64 _Mask);


__int64 _InterlockedOr64_np(__int64 volatile * _Value, __int64 _Mask);

char _InterlockedOr8(char volatile * _Value, char _Mask);


char _InterlockedOr8_np(char volatile * _Value, char _Mask);



long _InterlockedOr_np(long volatile * _Value, long _Mask);

long _InterlockedXor(long volatile * _Value, long _Mask);
short _InterlockedXor16(short volatile * _Value, short _Mask);


short _InterlockedXor16_np(short volatile * _Value, short _Mask);

__int64 _InterlockedXor64(__int64 volatile * _Value, __int64 _Mask);


__int64 _InterlockedXor64_np(__int64 volatile * _Value, __int64 _Mask);

char _InterlockedXor8(char volatile * _Value, char _Mask);


char _InterlockedXor8_np(char volatile * _Value, char _Mask);



long _InterlockedXor_np(long volatile * _Value, long _Mask);









void _ReadBarrier(void);






void _ReadWriteBarrier(void);
void * _ReturnAddress(void);

void _WriteBarrier(void);









void __addgsbyte(unsigned long, unsigned char);
void __addgsdword(unsigned long, unsigned long);
void __addgsqword(unsigned long, unsigned __int64);
void __addgsword(unsigned long, unsigned short);




void __code_seg(const char *);
void __cpuid(int[4], int);
void __cpuidex(int[4], int, int);
void __cdecl __debugbreak(void);

__int64 __emul(int, int);
unsigned __int64 __emulu(unsigned int, unsigned int);
__declspec(noreturn) void __fastfail(unsigned int);
void __faststorefence(void);
unsigned int __getcallerseflags(void);
void __halt(void);


unsigned char __inbyte(unsigned short);
void __inbytestring(unsigned short, unsigned char *, unsigned long);



void __incgsbyte(unsigned long);
void __incgsdword(unsigned long);
void __incgsqword(unsigned long);
void __incgsword(unsigned long);




unsigned long __indword(unsigned short);
void __indwordstring(unsigned short, unsigned long *, unsigned long);
void __int2c(void);
void __invlpg(void *);
unsigned short __inword(unsigned short);
void __inwordstring(unsigned short, unsigned short *, unsigned long);









void __lidt(void *);
unsigned __int64 __ll_lshift(unsigned __int64, int);
__int64 __ll_rshift(__int64, int);
unsigned int __lzcnt(unsigned int);
unsigned short __lzcnt16(unsigned short);
unsigned __int64 __lzcnt64(unsigned __int64);
void __movsb(unsigned char *, unsigned char const *, size_t);
void __movsd(unsigned long *, unsigned long const *, size_t);
void __movsq(unsigned long long *, unsigned long long const *, size_t);
void __movsw(unsigned short *, unsigned short const *, size_t);
__int64 __mulh(__int64, __int64);
void __nop(void);
void __nvreg_restore_fence(void);
void __nvreg_save_fence(void);
void __outbyte(unsigned short, unsigned char);
void __outbytestring(unsigned short, unsigned char *, unsigned long);
void __outdword(unsigned short, unsigned long);
void __outdwordstring(unsigned short, unsigned long *, unsigned long);
void __outword(unsigned short, unsigned short);
void __outwordstring(unsigned short, unsigned short *, unsigned long);
unsigned int __popcnt(unsigned int);
unsigned short __popcnt16(unsigned short);
unsigned __int64 __popcnt64(unsigned __int64);



unsigned __int64 __rdtsc(void);
unsigned __int64 __rdtscp(unsigned int *);
unsigned __int64 __readcr0(void);

unsigned __int64 __readcr2(void);

unsigned __int64 __readcr3(void);

unsigned __int64 __readcr4(void);

unsigned __int64 __readcr8(void);

unsigned __int64 __readdr(unsigned int);

unsigned __int64 __readeflags(void);





unsigned char __readgsbyte(unsigned long);
unsigned long __readgsdword(unsigned long);
unsigned __int64 __readgsqword(unsigned long);
unsigned short __readgsword(unsigned long);
unsigned __int64 __readmsr(unsigned long);
unsigned __int64 __readpmc(unsigned long);




unsigned long __segmentlimit(unsigned long);

unsigned __int64 __shiftleft128(unsigned __int64 _LowPart, unsigned __int64 _HighPart, unsigned char _Shift);
unsigned __int64 __shiftright128(unsigned __int64 _LowPart, unsigned __int64 _HighPart, unsigned char _Shift);
void __sidt(void *);

void __stosb(unsigned char *, unsigned char, size_t);
void __stosd(unsigned long *, unsigned long, size_t);
void __stosq(unsigned __int64 *, unsigned __int64, size_t);
void __stosw(unsigned short *, unsigned short, size_t);
void __svm_clgi(void);
void __svm_invlpga(void *, int);
void __svm_skinit(int);
void __svm_stgi(void);
void __svm_vmload(size_t);
void __svm_vmrun(size_t);
void __svm_vmsave(size_t);





void __ud2(void);
unsigned __int64 __ull_rshift(unsigned __int64, int);
unsigned __int64 __umulh(unsigned __int64, unsigned __int64);
void __vmx_off(void);
unsigned char __vmx_on(unsigned __int64 *);
unsigned char __vmx_vmclear(unsigned __int64 *);
unsigned char __vmx_vmlaunch(void);
unsigned char __vmx_vmptrld(unsigned __int64 *);
void __vmx_vmptrst(unsigned __int64 *);
unsigned char __vmx_vmread(size_t, size_t *);
unsigned char __vmx_vmresume(void);
unsigned char __vmx_vmwrite(size_t, size_t);
void __wbinvd(void);


void __writecr0(unsigned __int64);

void __writecr3(unsigned __int64);

void __writecr4(unsigned __int64);

void __writecr8(unsigned __int64);

void __writedr(unsigned int, unsigned __int64);

void __writeeflags(unsigned __int64);





void __writegsbyte(unsigned long, unsigned char);
void __writegsdword(unsigned long, unsigned long);
void __writegsqword(unsigned long, unsigned __int64);
void __writegsword(unsigned long, unsigned short);
void __writemsr(unsigned long, unsigned __int64);





unsigned char _bittest(long const *, long);
unsigned char _bittest64(__int64 const *, __int64);
unsigned char _bittestandcomplement(long *, long);
unsigned char _bittestandcomplement64(__int64 *, __int64);
unsigned char _bittestandreset(long *, long);
unsigned char _bittestandreset64(__int64 *, __int64);
unsigned char _bittestandset(long *, long);
unsigned char _bittestandset64(__int64 *, __int64);
  unsigned __int64 __cdecl _byteswap_uint64(  unsigned __int64);
  unsigned long __cdecl _byteswap_ulong(  unsigned long);
  unsigned short __cdecl _byteswap_ushort(  unsigned short);
void __cdecl _disable(void);
void __cdecl _enable(void);
unsigned char _interlockedbittestandreset(long volatile *, long);
unsigned char _interlockedbittestandreset64(__int64 volatile *, __int64);






unsigned char _interlockedbittestandset(long volatile *, long);
unsigned char _interlockedbittestandset64(__int64 volatile *, __int64);
















  unsigned long __cdecl _lrotl(  unsigned long,   int);
  unsigned long __cdecl _lrotr(  unsigned long,   int);




























































void _m_prefetch(void *);
void _m_prefetchw(volatile const void *);



































__m128i _mm_abs_epi16(__m128i);
__m128i _mm_abs_epi32(__m128i);
__m128i _mm_abs_epi8(__m128i);



__m128i _mm_add_epi16(__m128i, __m128i);
__m128i _mm_add_epi32(__m128i, __m128i);
__m128i _mm_add_epi64(__m128i, __m128i);
__m128i _mm_add_epi8(__m128i, __m128i);
__m128d _mm_add_pd(__m128d, __m128d);
__m128 _mm_add_ps(__m128, __m128);
__m128d _mm_add_sd(__m128d, __m128d);

__m128 _mm_add_ss(__m128, __m128);
__m128i _mm_adds_epi16(__m128i, __m128i);
__m128i _mm_adds_epi8(__m128i, __m128i);
__m128i _mm_adds_epu16(__m128i, __m128i);
__m128i _mm_adds_epu8(__m128i, __m128i);
__m128d _mm_addsub_pd(__m128d, __m128d);
__m128 _mm_addsub_ps(__m128, __m128);
__m128i _mm_alignr_epi8(__m128i, __m128i, int);

__m128d _mm_and_pd(__m128d, __m128d);
__m128 _mm_and_ps(__m128, __m128);
__m128i _mm_and_si128(__m128i, __m128i);
__m128d _mm_andnot_pd(__m128d, __m128d);
__m128 _mm_andnot_ps(__m128, __m128);
__m128i _mm_andnot_si128(__m128i, __m128i);
__m128i _mm_avg_epu16(__m128i, __m128i);
__m128i _mm_avg_epu8(__m128i, __m128i);
__m128i _mm_blend_epi16(__m128i, __m128i, int);
__m128d _mm_blend_pd(__m128d, __m128d, int);
__m128 _mm_blend_ps(__m128, __m128, int);
__m128i _mm_blendv_epi8(__m128i, __m128i, __m128i);
__m128d _mm_blendv_pd(__m128d, __m128d, __m128d);
__m128 _mm_blendv_ps(__m128, __m128, __m128);
void _mm_clflush(void const *);
void _mm_clflushopt(void const *);
void _mm_clwb(void const *);
void _mm_clzero(void const *);
__m128i _mm_cmpeq_epi16(__m128i, __m128i);
__m128i _mm_cmpeq_epi32(__m128i, __m128i);
__m128i _mm_cmpeq_epi64(__m128i, __m128i);
__m128i _mm_cmpeq_epi8(__m128i, __m128i);
__m128d _mm_cmpeq_pd(__m128d, __m128d);
__m128 _mm_cmpeq_ps(__m128, __m128);
__m128d _mm_cmpeq_sd(__m128d, __m128d);
__m128 _mm_cmpeq_ss(__m128, __m128);
int _mm_cmpestra(__m128i, int, __m128i, int, int);
int _mm_cmpestrc(__m128i, int, __m128i, int, int);
int _mm_cmpestri(__m128i, int, __m128i, int, int);
__m128i _mm_cmpestrm(__m128i, int, __m128i, int, int);
int _mm_cmpestro(__m128i, int, __m128i, int, int);
int _mm_cmpestrs(__m128i, int, __m128i, int, int);
int _mm_cmpestrz(__m128i, int, __m128i, int, int);
__m128d _mm_cmpge_pd(__m128d, __m128d);
__m128 _mm_cmpge_ps(__m128, __m128);
__m128d _mm_cmpge_sd(__m128d, __m128d);
__m128 _mm_cmpge_ss(__m128, __m128);
__m128i _mm_cmpgt_epi16(__m128i, __m128i);
__m128i _mm_cmpgt_epi32(__m128i, __m128i);
__m128i _mm_cmpgt_epi64(__m128i, __m128i);
__m128i _mm_cmpgt_epi8(__m128i, __m128i);
__m128d _mm_cmpgt_pd(__m128d, __m128d);
__m128 _mm_cmpgt_ps(__m128, __m128);
__m128d _mm_cmpgt_sd(__m128d, __m128d);
__m128 _mm_cmpgt_ss(__m128, __m128);
int _mm_cmpistra(__m128i, __m128i, int);
int _mm_cmpistrc(__m128i, __m128i, int);
int _mm_cmpistri(__m128i, __m128i, int);
__m128i _mm_cmpistrm(__m128i, __m128i, int);
int _mm_cmpistro(__m128i, __m128i, int);
int _mm_cmpistrs(__m128i, __m128i, int);
int _mm_cmpistrz(__m128i, __m128i, int);
__m128d _mm_cmple_pd(__m128d, __m128d);
__m128 _mm_cmple_ps(__m128, __m128);
__m128d _mm_cmple_sd(__m128d, __m128d);
__m128 _mm_cmple_ss(__m128, __m128);
__m128i _mm_cmplt_epi16(__m128i, __m128i);
__m128i _mm_cmplt_epi32(__m128i, __m128i);
__m128i _mm_cmplt_epi8(__m128i, __m128i);
__m128d _mm_cmplt_pd(__m128d, __m128d);
__m128 _mm_cmplt_ps(__m128, __m128);
__m128d _mm_cmplt_sd(__m128d, __m128d);
__m128 _mm_cmplt_ss(__m128, __m128);
__m128d _mm_cmpneq_pd(__m128d, __m128d);
__m128 _mm_cmpneq_ps(__m128, __m128);
__m128d _mm_cmpneq_sd(__m128d, __m128d);
__m128 _mm_cmpneq_ss(__m128, __m128);
__m128d _mm_cmpnge_pd(__m128d, __m128d);
__m128 _mm_cmpnge_ps(__m128, __m128);
__m128d _mm_cmpnge_sd(__m128d, __m128d);
__m128 _mm_cmpnge_ss(__m128, __m128);
__m128d _mm_cmpngt_pd(__m128d, __m128d);
__m128 _mm_cmpngt_ps(__m128, __m128);
__m128d _mm_cmpngt_sd(__m128d, __m128d);
__m128 _mm_cmpngt_ss(__m128, __m128);
__m128d _mm_cmpnle_pd(__m128d, __m128d);
__m128 _mm_cmpnle_ps(__m128, __m128);
__m128d _mm_cmpnle_sd(__m128d, __m128d);
__m128 _mm_cmpnle_ss(__m128, __m128);
__m128d _mm_cmpnlt_pd(__m128d, __m128d);
__m128 _mm_cmpnlt_ps(__m128, __m128);
__m128d _mm_cmpnlt_sd(__m128d, __m128d);
__m128 _mm_cmpnlt_ss(__m128, __m128);
__m128d _mm_cmpord_pd(__m128d, __m128d);
__m128 _mm_cmpord_ps(__m128, __m128);
__m128d _mm_cmpord_sd(__m128d, __m128d);
__m128 _mm_cmpord_ss(__m128, __m128);
__m128d _mm_cmpunord_pd(__m128d, __m128d);
__m128 _mm_cmpunord_ps(__m128, __m128);
__m128d _mm_cmpunord_sd(__m128d, __m128d);
__m128 _mm_cmpunord_ss(__m128, __m128);
int _mm_comieq_sd(__m128d, __m128d);
int _mm_comieq_ss(__m128, __m128);
int _mm_comige_sd(__m128d, __m128d);
int _mm_comige_ss(__m128, __m128);
int _mm_comigt_sd(__m128d, __m128d);
int _mm_comigt_ss(__m128, __m128);
int _mm_comile_sd(__m128d, __m128d);
int _mm_comile_ss(__m128, __m128);
int _mm_comilt_sd(__m128d, __m128d);
int _mm_comilt_ss(__m128, __m128);
int _mm_comineq_sd(__m128d, __m128d);
int _mm_comineq_ss(__m128, __m128);
unsigned int _mm_crc32_u16(unsigned int, unsigned short);
unsigned int _mm_crc32_u32(unsigned int, unsigned int);
unsigned __int64 _mm_crc32_u64(unsigned __int64, unsigned __int64);
unsigned int _mm_crc32_u8(unsigned int, unsigned char);


__m128 _mm_cvt_si2ss(__m128, int);
int _mm_cvt_ss2si(__m128);
__m128i _mm_cvtepi16_epi32(__m128i);
__m128i _mm_cvtepi16_epi64(__m128i);
__m128i _mm_cvtepi32_epi64(__m128i);
__m128d _mm_cvtepi32_pd(__m128i);
__m128 _mm_cvtepi32_ps(__m128i);
__m128i _mm_cvtepi8_epi16(__m128i);
__m128i _mm_cvtepi8_epi32(__m128i);
__m128i _mm_cvtepi8_epi64(__m128i);
__m128i _mm_cvtepu16_epi32(__m128i);
__m128i _mm_cvtepu16_epi64(__m128i);
__m128i _mm_cvtepu32_epi64(__m128i);
__m128i _mm_cvtepu8_epi16(__m128i);
__m128i _mm_cvtepu8_epi32(__m128i);
__m128i _mm_cvtepu8_epi64(__m128i);
__m128i _mm_cvtpd_epi32(__m128d);

__m128 _mm_cvtpd_ps(__m128d);

__m128i _mm_cvtps_epi32(__m128);
__m128d _mm_cvtps_pd(__m128);
int _mm_cvtsd_si32(__m128d);
__int64 _mm_cvtsd_si64(__m128d);
__int64 _mm_cvtsd_si64x(__m128d);
__m128 _mm_cvtsd_ss(__m128, __m128d);
int _mm_cvtsi128_si32(__m128i);
__int64 _mm_cvtsi128_si64(__m128i);
__int64 _mm_cvtsi128_si64x(__m128i);
__m128d _mm_cvtsi32_sd(__m128d, int);
__m128i _mm_cvtsi32_si128(int);
__m128d _mm_cvtsi64_sd(__m128d, __int64);
__m128i _mm_cvtsi64_si128(__int64);
__m128 _mm_cvtsi64_ss(__m128, __int64);
__m128d _mm_cvtsi64x_sd(__m128d, __int64);
__m128i _mm_cvtsi64x_si128(__int64);
__m128 _mm_cvtsi64x_ss(__m128, __int64);
__m128d _mm_cvtss_sd(__m128d, __m128);
__int64 _mm_cvtss_si64(__m128);
__int64 _mm_cvtss_si64x(__m128);

int _mm_cvtt_ss2si(__m128);
__m128i _mm_cvttpd_epi32(__m128d);

__m128i _mm_cvttps_epi32(__m128);
int _mm_cvttsd_si32(__m128d);
__int64 _mm_cvttsd_si64(__m128d);
__int64 _mm_cvttsd_si64x(__m128d);
__int64 _mm_cvttss_si64(__m128);
__int64 _mm_cvttss_si64x(__m128);
__m128d _mm_div_pd(__m128d, __m128d);
__m128 _mm_div_ps(__m128, __m128);
__m128d _mm_div_sd(__m128d, __m128d);
__m128 _mm_div_ss(__m128, __m128);
__m128d _mm_dp_pd(__m128d, __m128d, int);
__m128 _mm_dp_ps(__m128, __m128, int);
int _mm_extract_epi16(__m128i, int);
int _mm_extract_epi32(__m128i, int);
__int64 _mm_extract_epi64(__m128i, int);
int _mm_extract_epi8(__m128i, int);
int _mm_extract_ps(__m128, int);
__m128i _mm_extract_si64(__m128i, __m128i);
__m128i _mm_extracti_si64(__m128i, int, int);
unsigned int _mm_getcsr(void);
__m128i _mm_hadd_epi16(__m128i, __m128i);
__m128i _mm_hadd_epi32(__m128i, __m128i);
__m128d _mm_hadd_pd(__m128d, __m128d);


__m128 _mm_hadd_ps(__m128, __m128);
__m128i _mm_hadds_epi16(__m128i, __m128i);

__m128i _mm_hsub_epi16(__m128i, __m128i);
__m128i _mm_hsub_epi32(__m128i, __m128i);
__m128d _mm_hsub_pd(__m128d, __m128d);


__m128 _mm_hsub_ps(__m128, __m128);
__m128i _mm_hsubs_epi16(__m128i, __m128i);

__m128i _mm_insert_epi16(__m128i, int, int);
__m128i _mm_insert_epi32(__m128i, int, int);
__m128i _mm_insert_epi64(__m128i, __int64, int);
__m128i _mm_insert_epi8(__m128i, int, int);
__m128 _mm_insert_ps(__m128, __m128, int);
__m128i _mm_insert_si64(__m128i, __m128i);
__m128i _mm_inserti_si64(__m128i, __m128i, int, int);
__m128i _mm_lddqu_si128(__m128i const *);
void _mm_lfence(void);
__m128d _mm_load1_pd(double const *);
__m128d _mm_load_pd(double const *);
__m128 _mm_load_ps(float const *);
__m128 _mm_load_ps1(float const *);
__m128d _mm_load_sd(double const *);
__m128i _mm_load_si128(__m128i const *);
__m128 _mm_load_ss(float const *);
__m128d _mm_loaddup_pd(double const *);
__m128d _mm_loadh_pd(__m128d, double const *);
__m128 _mm_loadh_pi(__m128, __m64 const *);
__m128i _mm_loadl_epi64(__m128i const *);
__m128d _mm_loadl_pd(__m128d, double const *);
__m128 _mm_loadl_pi(__m128, __m64 const *);
__m128d _mm_loadr_pd(double const *);
__m128 _mm_loadr_ps(float const *);
__m128d _mm_loadu_pd(double const *);
__m128 _mm_loadu_ps(float const *);
__m128i _mm_loadu_si128(__m128i const *);
__m128i _mm_madd_epi16(__m128i, __m128i);
__m128i _mm_maddubs_epi16(__m128i, __m128i);

void _mm_maskmoveu_si128(__m128i, __m128i, char *);
__m128i _mm_max_epi16(__m128i, __m128i);
__m128i _mm_max_epi32(__m128i, __m128i);
__m128i _mm_max_epi8(__m128i, __m128i);
__m128i _mm_max_epu16(__m128i, __m128i);
__m128i _mm_max_epu32(__m128i, __m128i);
__m128i _mm_max_epu8(__m128i, __m128i);
__m128d _mm_max_pd(__m128d, __m128d);
__m128 _mm_max_ps(__m128, __m128);
__m128d _mm_max_sd(__m128d, __m128d);
__m128 _mm_max_ss(__m128, __m128);
void _mm_mfence(void);
__m128i _mm_min_epi16(__m128i, __m128i);
__m128i _mm_min_epi32(__m128i, __m128i);
__m128i _mm_min_epi8(__m128i, __m128i);
__m128i _mm_min_epu16(__m128i, __m128i);
__m128i _mm_min_epu32(__m128i, __m128i);
__m128i _mm_min_epu8(__m128i, __m128i);
__m128d _mm_min_pd(__m128d, __m128d);
__m128 _mm_min_ps(__m128, __m128);
__m128d _mm_min_sd(__m128d, __m128d);
__m128 _mm_min_ss(__m128, __m128);
__m128i _mm_minpos_epu16(__m128i);
void _mm_monitor(void const *, unsigned int, unsigned int);
__m128i _mm_move_epi64(__m128i);
__m128d _mm_move_sd(__m128d, __m128d);
__m128 _mm_move_ss(__m128, __m128);
__m128d _mm_movedup_pd(__m128d);
__m128 _mm_movehdup_ps(__m128);
__m128 _mm_movehl_ps(__m128, __m128);
__m128 _mm_moveldup_ps(__m128);
__m128 _mm_movelh_ps(__m128, __m128);
int _mm_movemask_epi8(__m128i);
int _mm_movemask_pd(__m128d);
int _mm_movemask_ps(__m128);


__m128i _mm_mpsadbw_epu8(__m128i, __m128i, int);
__m128i _mm_mul_epi32(__m128i, __m128i);
__m128i _mm_mul_epu32(__m128i, __m128i);
__m128d _mm_mul_pd(__m128d, __m128d);
__m128 _mm_mul_ps(__m128, __m128);
__m128d _mm_mul_sd(__m128d, __m128d);
__m128 _mm_mul_ss(__m128, __m128);

__m128i _mm_mulhi_epi16(__m128i, __m128i);
__m128i _mm_mulhi_epu16(__m128i, __m128i);
__m128i _mm_mulhrs_epi16(__m128i, __m128i);

__m128i _mm_mullo_epi16(__m128i, __m128i);
__m128i _mm_mullo_epi32(__m128i, __m128i);
void _mm_mwait(unsigned int, unsigned int);
__m128d _mm_or_pd(__m128d, __m128d);
__m128 _mm_or_ps(__m128, __m128);
__m128i _mm_or_si128(__m128i, __m128i);
__m128i _mm_packs_epi16(__m128i, __m128i);
__m128i _mm_packs_epi32(__m128i, __m128i);
__m128i _mm_packus_epi16(__m128i, __m128i);
__m128i _mm_packus_epi32(__m128i, __m128i);
void _mm_pause(void);
int _mm_popcnt_u32(unsigned int);
__int64 _mm_popcnt_u64(unsigned __int64);
void _mm_prefetch(char const *, int);
__m128 _mm_rcp_ps(__m128);
__m128 _mm_rcp_ss(__m128);
__m128d _mm_round_pd(__m128d, int);
__m128 _mm_round_ps(__m128, int);
__m128d _mm_round_sd(__m128d, __m128d, int);
__m128 _mm_round_ss(__m128, __m128, int);
__m128 _mm_rsqrt_ps(__m128);
__m128 _mm_rsqrt_ss(__m128);
__m128i _mm_sad_epu8(__m128i, __m128i);
__m128i _mm_set1_epi16(short);
__m128i _mm_set1_epi32(int);

__m128i _mm_set1_epi64x(__int64);
__m128i _mm_set1_epi8(char);
__m128d _mm_set1_pd(double);



__m128i _mm_set_epi16(short, short, short, short, short, short, short, short);
__m128i _mm_set_epi32(int, int, int, int);

__m128i _mm_set_epi64x(__int64, __int64);
__m128i _mm_set_epi8(char, char, char, char, char, char, char, char, char, char, char, char, char, char, char, char);
__m128d _mm_set_pd(double, double);



__m128 _mm_set_ps(float, float, float, float);
__m128 _mm_set_ps1(float);
__m128d _mm_set_sd(double);
__m128 _mm_set_ss(float);
void _mm_setcsr(unsigned int);
__m128i _mm_setl_epi64(__m128i);
__m128i _mm_setr_epi16(short, short, short, short, short, short, short, short);
__m128i _mm_setr_epi32(int, int, int, int);

__m128i _mm_setr_epi64x(__int64, __int64);
__m128i _mm_setr_epi8(char, char, char, char, char, char, char, char, char, char, char, char, char, char, char, char);
__m128d _mm_setr_pd(double, double);



__m128 _mm_setr_ps(float, float, float, float);
__m128d _mm_setzero_pd(void);
__m128 _mm_setzero_ps(void);
__m128i _mm_setzero_si128(void);

void _mm_sfence(void);
__m128i _mm_shuffle_epi32(__m128i, int);
__m128i _mm_shuffle_epi8(__m128i, __m128i);
__m128d _mm_shuffle_pd(__m128d, __m128d, int);

__m128 _mm_shuffle_ps(__m128, __m128, unsigned int);
__m128i _mm_shufflehi_epi16(__m128i, int);
__m128i _mm_shufflelo_epi16(__m128i, int);
__m128i _mm_sign_epi16(__m128i, __m128i);
__m128i _mm_sign_epi32(__m128i, __m128i);
__m128i _mm_sign_epi8(__m128i, __m128i);



__m128i _mm_sll_epi16(__m128i, __m128i);
__m128i _mm_sll_epi32(__m128i, __m128i);
__m128i _mm_sll_epi64(__m128i, __m128i);
__m128i _mm_slli_epi16(__m128i, int);
__m128i _mm_slli_epi32(__m128i, int);
__m128i _mm_slli_epi64(__m128i, int);
__m128i _mm_slli_si128(__m128i, int);
__m128d _mm_sqrt_pd(__m128d);
__m128 _mm_sqrt_ps(__m128);
__m128d _mm_sqrt_sd(__m128d, __m128d);
__m128 _mm_sqrt_ss(__m128);
__m128i _mm_sra_epi16(__m128i, __m128i);
__m128i _mm_sra_epi32(__m128i, __m128i);
__m128i _mm_srai_epi16(__m128i, int);
__m128i _mm_srai_epi32(__m128i, int);
__m128i _mm_srl_epi16(__m128i, __m128i);
__m128i _mm_srl_epi32(__m128i, __m128i);
__m128i _mm_srl_epi64(__m128i, __m128i);
__m128i _mm_srli_epi16(__m128i, int);
__m128i _mm_srli_epi32(__m128i, int);
__m128i _mm_srli_epi64(__m128i, int);
__m128i _mm_srli_si128(__m128i, int);
void _mm_store1_pd(double *, __m128d);
void _mm_store_pd(double *, __m128d);
void _mm_store_ps(float *, __m128);
void _mm_store_ps1(float *, __m128);
void _mm_store_sd(double *, __m128d);
void _mm_store_si128(__m128i *, __m128i);
void _mm_store_ss(float *, __m128);
void _mm_storeh_pd(double *, __m128d);
void _mm_storeh_pi(__m64 *, __m128);
void _mm_storel_epi64(__m128i *, __m128i);
void _mm_storel_pd(double *, __m128d);
void _mm_storel_pi(__m64 *, __m128);
void _mm_storer_pd(double *, __m128d);
void _mm_storer_ps(float *, __m128);
void _mm_storeu_pd(double *, __m128d);
void _mm_storeu_ps(float *, __m128);
void _mm_storeu_si128(__m128i *, __m128i);
__m128i _mm_stream_load_si128(const __m128i *);
void _mm_stream_pd(double *, __m128d);

void _mm_stream_ps(float *, __m128);
void _mm_stream_sd(double *, __m128d);
void _mm_stream_si128(__m128i *, __m128i);
void _mm_stream_si32(int *, int);
void _mm_stream_si64x(__int64 *, __int64);
void _mm_stream_ss(float *, __m128);
__m128i _mm_sub_epi16(__m128i, __m128i);
__m128i _mm_sub_epi32(__m128i, __m128i);
__m128i _mm_sub_epi64(__m128i, __m128i);
__m128i _mm_sub_epi8(__m128i, __m128i);
__m128d _mm_sub_pd(__m128d, __m128d);
__m128 _mm_sub_ps(__m128, __m128);
__m128d _mm_sub_sd(__m128d, __m128d);

__m128 _mm_sub_ss(__m128, __m128);
__m128i _mm_subs_epi16(__m128i, __m128i);
__m128i _mm_subs_epi8(__m128i, __m128i);
__m128i _mm_subs_epu16(__m128i, __m128i);
__m128i _mm_subs_epu8(__m128i, __m128i);
int _mm_testc_si128(__m128i, __m128i);
int _mm_testnzc_si128(__m128i, __m128i);
int _mm_testz_si128(__m128i, __m128i);
int _mm_ucomieq_sd(__m128d, __m128d);
int _mm_ucomieq_ss(__m128, __m128);
int _mm_ucomige_sd(__m128d, __m128d);
int _mm_ucomige_ss(__m128, __m128);
int _mm_ucomigt_sd(__m128d, __m128d);
int _mm_ucomigt_ss(__m128, __m128);
int _mm_ucomile_sd(__m128d, __m128d);
int _mm_ucomile_ss(__m128, __m128);
int _mm_ucomilt_sd(__m128d, __m128d);
int _mm_ucomilt_ss(__m128, __m128);
int _mm_ucomineq_sd(__m128d, __m128d);
int _mm_ucomineq_ss(__m128, __m128);
__m128i _mm_unpackhi_epi16(__m128i, __m128i);
__m128i _mm_unpackhi_epi32(__m128i, __m128i);
__m128i _mm_unpackhi_epi64(__m128i, __m128i);
__m128i _mm_unpackhi_epi8(__m128i, __m128i);
__m128d _mm_unpackhi_pd(__m128d, __m128d);
__m128 _mm_unpackhi_ps(__m128, __m128);
__m128i _mm_unpacklo_epi16(__m128i, __m128i);
__m128i _mm_unpacklo_epi32(__m128i, __m128i);
__m128i _mm_unpacklo_epi64(__m128i, __m128i);
__m128i _mm_unpacklo_epi8(__m128i, __m128i);
__m128d _mm_unpacklo_pd(__m128d, __m128d);
__m128 _mm_unpacklo_ps(__m128, __m128);
__m128d _mm_xor_pd(__m128d, __m128d);
__m128 _mm_xor_ps(__m128, __m128);
__m128i _mm_xor_si128(__m128i, __m128i);
__int64 _mul128(__int64 _Multiplier, __int64 _Multiplicand, __int64 * _HighProduct);
unsigned int __cdecl _rotl(  unsigned int _Value,   int _Shift);
unsigned short __cdecl _rotl16(unsigned short _Value, unsigned char _Shift);
unsigned __int64 __cdecl _rotl64(  unsigned __int64 _Value,   int _Shift);
unsigned char __cdecl _rotl8(unsigned char _Value, unsigned char _Shift);
unsigned int __cdecl _rotr(  unsigned int _Value,   int _Shift);
unsigned short __cdecl _rotr16(unsigned short _Value, unsigned char _Shift);
unsigned __int64 __cdecl _rotr64(  unsigned __int64 _Value,   int _Shift);
unsigned char __cdecl _rotr8(unsigned char _Value, unsigned char _Shift);
int __cdecl _setjmp(jmp_buf);
int __cdecl _setjmpex(jmp_buf);
unsigned __int64 _umul128(unsigned __int64 _Multiplier, unsigned __int64 _Multiplicand, unsigned __int64 * _HighProduct);
void _rsm(void);
void _lgdt(void *);
void _sgdt(void *);
void _clac(void);
void _stac(void);
unsigned char __cdecl _addcarry_u8(unsigned char, unsigned char, unsigned char, unsigned char *);
unsigned char __cdecl _subborrow_u8(unsigned char, unsigned char, unsigned char, unsigned char *);
unsigned char __cdecl _addcarry_u16(unsigned char, unsigned short, unsigned short, unsigned short *);
unsigned char __cdecl _subborrow_u16(unsigned char, unsigned short, unsigned short, unsigned short *);
unsigned char __cdecl _addcarry_u32(unsigned char, unsigned int, unsigned int, unsigned int *);
unsigned char __cdecl _subborrow_u32(unsigned char, unsigned int, unsigned int, unsigned int *);
unsigned char __cdecl _addcarry_u64(unsigned char, unsigned __int64, unsigned __int64, unsigned __int64 *);
unsigned char __cdecl _subborrow_u64(unsigned char, unsigned __int64, unsigned __int64, unsigned __int64 *);
void _mm_monitorx(void const *, unsigned int, unsigned int);
void _mm_mwaitx(unsigned int, unsigned int, unsigned int);


}






   

   


  








 
 #pragma warning(pop)
 #pragma pack(pop)










 #pragma pack(push,8)
 #pragma warning(push,3)
 
 

 #pragma warning(disable: 4700)

namespace std {
		
template<class _Ty> inline
	pair<_Ty *, ptrdiff_t>
		get_temporary_buffer(ptrdiff_t _Count) noexcept
	{	
	_Ty *_Pbuf;

	if (_Count < 0)
		_Count = 0;
	else if (((size_t)(-1) / sizeof (_Ty) < _Count))
		_Xbad_alloc();	
	for (_Pbuf = 0; 0 < _Count; _Count /= 2)
		if ((_Pbuf = (_Ty *)operator new(
			(size_t)_Count * sizeof (_Ty), nothrow)) != 0)
			break;

	return (pair<_Ty *, ptrdiff_t>(_Pbuf, _Count));
	}

		
template<class _Ty> inline
	void return_temporary_buffer(_Ty *_Pbuf)
	{	
	operator delete(_Pbuf);
	}

		
template<class _InIt,
	class _FwdIt> inline
	_FwdIt _Uninitialized_copy_unchecked1(_InIt _First, _InIt _Last,
		_FwdIt _Dest, _General_ptr_iterator_tag)
	{	
	_FwdIt _Next = _Dest;

	try {
	for (; _First != _Last; ++_Dest, (void)++_First)
		_Construct(_Unfancy(_Dest), *_First);
	} catch (...) {
	_Destroy_range(_Next, _Dest);
	throw;
	}

	return (_Dest);
	}

template<class _InIt,
	class _FwdIt> inline
	_FwdIt _Uninitialized_copy_unchecked1(_InIt _First, _InIt _Last,
		_FwdIt _Dest, _Really_trivial_ptr_iterator_tag)
	{	
	return (_Copy_memmove(_First, _Last, _Dest));
	}

template<class _InIt,
	class _FwdIt> inline
	_FwdIt _Uninitialized_copy_unchecked(_InIt _First, _InIt _Last,
		_FwdIt _Dest)
	{	
	return (_Uninitialized_copy_unchecked1(_First, _Last,
		_Dest, _Ptr_copy_cat(_First, _Dest)));
	}

template<class _InIt,
	class _FwdIt> inline
	_FwdIt _Uninitialized_copy1(_InIt _First, _InIt _Last,
		_FwdIt _Dest, input_iterator_tag, forward_iterator_tag)
	{	
	return (_Rechecked(_Dest,
		_Uninitialized_copy_unchecked(_First, _Last, _Unchecked_idl0(_Dest))));
	}

template<class _InIt,
	class _FwdIt> inline
	_FwdIt _Uninitialized_copy1(_InIt _First, _InIt _Last,
		_FwdIt _Dest, random_access_iterator_tag, random_access_iterator_tag)
	{	
	;
	return (_Rechecked(_Dest,
		_Uninitialized_copy_unchecked(_First, _Last, _Unchecked(_Dest))));
	}

template<class _InIt,
	class _FwdIt> inline
	_FwdIt uninitialized_copy(_InIt _First, _InIt _Last,
		_FwdIt _Dest)
	{	
	struct _Unchecked_iterators { static void  _Deprecate(false_type) { } static void _Deprecate(true_type) { } }; (_Unchecked_iterators::_Deprecate(_Is_checked(_Dest)));
	;
	return (_Uninitialized_copy1(_Unchecked(_First), _Unchecked(_Last),
		_Dest, _Iter_cat_t<_InIt>(), _Iter_cat_t<_FwdIt>()));
	}

 












		
template<class _InIt,
	class _Diff,
	class _FwdIt> inline
	_FwdIt _Uninitialized_copy_n_unchecked1(_InIt _First, _Diff _Count,
		_FwdIt _Dest, _General_ptr_iterator_tag)
	{	
	_FwdIt _Next = _Dest;

	try {
	for (; 0 < _Count; --_Count, (void)++_Dest, ++_First)
		_Construct(_Unfancy(_Dest), *_First);
	} catch (...) {
	_Destroy_range(_Next, _Dest);
	throw;
	}

	return (_Dest);
	}

template<class _InIt,
	class _Diff,
	class _FwdIt> inline
	_FwdIt _Uninitialized_copy_n_unchecked1(_InIt _First, _Diff _Count,
		_FwdIt _Dest, _Really_trivial_ptr_iterator_tag)
	{	
	if (0 < _Count)
		return (_Copy_memmove(_First, _First + _Count, _Dest));
	return (_Dest);
	}

template<class _InIt,
	class _Diff,
	class _FwdIt> inline
	_FwdIt _Uninitialized_copy_n_unchecked(_InIt _First, _Diff _Count,
		_FwdIt _Dest)
	{	
	return (_Uninitialized_copy_n_unchecked1(_First, _Count,
		_Dest, _Ptr_copy_cat(_First, _Dest)));
	}

template<class _InIt,
	class _Diff,
	class _FwdIt> inline
	_FwdIt uninitialized_copy_n(_InIt _First, _Diff _Count,
		_FwdIt _Dest)
	{	
		
	struct _Unchecked_iterators { static void  _Deprecate(false_type) { } static void _Deprecate(true_type) { } }; (_Unchecked_iterators::_Deprecate(_Is_checked(_Dest)));
	return (_Rechecked(_Dest,
		_Uninitialized_copy_n_unchecked(_Unchecked_n(_First, _Count), _Count, _Unchecked_n(_Dest, _Count))));
	}

 







































		
template<class _InIt,
	class _FwdIt,
	class _Alloc> inline
	_FwdIt _Uninitialized_copy_al_unchecked1(_InIt _First, _InIt _Last, _FwdIt _Dest,
		_Wrap_alloc<_Alloc>& _Al, _General_ptr_iterator_tag, _Any_tag)
	{	
	_FwdIt _Next = _Dest;

	try {
	for (; _First != _Last; ++_Dest, (void)++_First)
		_Al.construct(_Unfancy(_Dest), *_First);
	} catch (...) {
	_Destroy_range(_Next, _Dest, _Al);
	throw;
	}

	return (_Dest);
	}

template<class _Ty1,
	class _Ty2,
	class _Alloc> inline
	_Ty2 *_Uninitialized_copy_al_unchecked1(_Ty1 *_First, _Ty1 *_Last, _Ty2 *_Dest,
		_Wrap_alloc<_Alloc>&, _Really_trivial_ptr_iterator_tag, true_type)
	{	
	return (_Copy_memmove(_First, _Last, _Dest));
	}

template<class _InIt,
	class _FwdIt,
	class _Alloc> inline
	_FwdIt _Uninitialized_copy_al_unchecked(_InIt _First, _InIt _Last, _FwdIt _Dest,
		_Wrap_alloc<_Alloc>& _Al)
	{	
	return (_Uninitialized_copy_al_unchecked1(_First, _Last, _Dest, _Al,
		_Ptr_copy_cat(_First, _Dest),
		_Uses_default_construct_t<_Alloc, decltype(_Unfancy(_Dest)), decltype(*_First)>()));
	}

template<class _InIt,
	class _FwdIt,
	class _Alloc> inline
	_FwdIt _Uninitialized_copy(_InIt _First, _InIt _Last, _FwdIt _Dest,
		_Wrap_alloc<_Alloc>& _Al)
	{	
		
		
	return (_Rechecked(_Dest,
		_Uninitialized_copy_al_unchecked(_Unchecked(_First), _Unchecked(_Last),
			_Unchecked(_Dest), _Al)));
	}

		
template<class _InIt,
	class _FwdIt,
	class _Alloc> inline
	_FwdIt _Uninitialized_move_al_unchecked1(_InIt _First, _InIt _Last, _FwdIt _Dest,
		_Wrap_alloc<_Alloc>& _Al, _General_ptr_iterator_tag, _Any_tag)
	{	
	_FwdIt _Next = _Dest;

	try {
	for (; _First != _Last; ++_Dest, (void)++_First)
		_Al.construct(_Unfancy(_Dest), ::std:: move(*_First));
	} catch (...) {
	_Destroy_range(_Next, _Dest, _Al);
	throw;
	}

	return (_Dest);
	}

template<class _Ty1,
	class _Ty2,
	class _Alloc> inline
	_Ty2 *_Uninitialized_move_al_unchecked1(_Ty1 *_First, _Ty1 *_Last, _Ty2 *_Dest,
		_Wrap_alloc<_Alloc>&, _Really_trivial_ptr_iterator_tag, true_type)
	{	
	return (_Copy_memmove(_First, _Last, _Dest));
	}

template<class _InIt,
	class _FwdIt,
	class _Alloc> inline
	_FwdIt _Uninitialized_move_al_unchecked(_InIt _First, _InIt _Last, _FwdIt _Dest,
		_Wrap_alloc<_Alloc>& _Al)
	{	
	typedef decltype(::std:: move(*_First)) _Src_type; 
	return (_Uninitialized_move_al_unchecked1(_First, _Last, _Dest, _Al,
		_Ptr_move_cat(_First, _Dest),
		_Uses_default_construct_t<_Alloc, decltype(_Unfancy(_Dest)), _Src_type>()));
	}

template<class _InIt,
	class _FwdIt,
	class _Alloc> inline
	_FwdIt _Uninitialized_move(_InIt _First, _InIt _Last, _FwdIt _Dest,
		_Wrap_alloc<_Alloc>& _Al)
	{	
		
		
	return (_Rechecked(_Dest,
		_Uninitialized_move_al_unchecked(_Unchecked(_First), _Unchecked(_Last),
			_Unchecked(_Dest), _Al)));
	}

		
template<class _FwdIt,
	class _Tval> inline
	void _Uninitialized_fill_unchecked1(_FwdIt _First, _FwdIt _Last, const _Tval& _Val, false_type)
	{	
	_FwdIt _Next = _First;

	try {
	for (; _First != _Last; ++_First)
		_Construct(_Unfancy(_First), _Val);
	} catch (...) {
	_Destroy_range(_Next, _First);
	throw;
	}
	}

template<class _FwdIt,
	class _Tval> inline
	void _Uninitialized_fill_unchecked1(_FwdIt _First, _FwdIt _Last, const _Tval& _Val, true_type)
	{	
	:: memset(_First, _Val, _Last - _First);
	}

template<class _FwdIt,
	class _Tval> inline
	void _Uninitialized_fill_unchecked(_FwdIt _First, _FwdIt _Last, const _Tval& _Val)
	{	
	_Uninitialized_fill_unchecked1(_First, _Last, _Val, _Fill_memset_is_safe(_First, _Val));
	}

template<class _FwdIt,
	class _Tval> inline
	void uninitialized_fill(_FwdIt _First, _FwdIt _Last, const _Tval& _Val)
	{	
	;
	_Uninitialized_fill_unchecked(_Unchecked(_First), _Unchecked(_Last), _Val);
	}

		
template<class _FwdIt,
	class _Diff,
	class _Tval> inline
	_FwdIt _Uninitialized_fill_n_unchecked1(_FwdIt _First, _Diff _Count, const _Tval& _Val, false_type)
	{	
	_FwdIt _Next = _First;

	try {
	for (; 0 < _Count; --_Count, (void)++_First)
		_Construct(_Unfancy(_First), _Val);
	} catch (...) {
	_Destroy_range(_Next, _First);
	throw;
	}

	return (_First);
	}

template<class _FwdIt,
	class _Diff,
	class _Tval> inline
	_FwdIt _Uninitialized_fill_n_unchecked1(_FwdIt _First, _Diff _Count, const _Tval& _Val, true_type)
	{	
	if (0 < _Count)
		{
		:: memset(_First, _Val, _Count);
		return (_First + _Count);
		}

	return (_First);
	}

template<class _FwdIt,
	class _Diff,
	class _Tval> inline
	_FwdIt _Uninitialized_fill_n_unchecked(_FwdIt _First, _Diff _Count, const _Tval& _Val)
	{	
	return (_Uninitialized_fill_n_unchecked1(_First, _Count, _Val, _Fill_memset_is_safe(_First, _Val)));
	}

template<class _FwdIt,
	class _Diff,
	class _Tval> inline
	_FwdIt uninitialized_fill_n(_FwdIt _First, _Diff _Count,
		const _Tval& _Val)
	{	
	return (_Rechecked(_First,
		_Uninitialized_fill_n_unchecked(_Unchecked_n(_First, _Count), _Count, _Val)));
	}

		
template<class _FwdIt,
	class _Diff,
	class _Alloc> inline
	void _Uninit_alloc_fill_n1(_FwdIt _First, _Diff _Count, const _Iter_value_t<_FwdIt> * _Pval,
		_Wrap_alloc<_Alloc>& _Al, false_type)
	{	
	_FwdIt _Next = _First;

	try {
	for (; 0 < _Count; --_Count, (void)++_First)
		_Al.construct(_Unfancy(_First), *_Pval);
	} catch (...) {
	_Destroy_range(_Next, _First, _Al);
	throw;
	}
	}

template<class _FwdIt,
	class _Diff,
	class _Alloc> inline
	void _Uninit_alloc_fill_n1(_FwdIt _First, _Diff _Count, const _Iter_value_t<_FwdIt> * _Pval,
		_Wrap_alloc<_Alloc>&, true_type)
	{	
	:: memset(_First, *_Pval, _Count);
	}

template<class _FwdIt,
	class _Diff,
	class _Alloc> inline
	void _Uninitialized_fill_n(_FwdIt _First, _Diff _Count,
		const _Iter_value_t<_FwdIt> * _Pval, _Wrap_alloc<_Alloc>& _Al)
	{	
	_Uninit_alloc_fill_n1(_First, _Count, _Pval, _Al,
		typename conjunction<decltype(_Fill_memset_is_safe(_First, *_Pval)),
			_Uses_default_construct<_Alloc, decltype(_Unfancy(_First)), decltype(*_Pval)>>::type());
	}

template<class _FwdIt,
	class _Diff,
	class _Alloc> inline
	void _Uninitialized_default_fill_n1(_FwdIt _First, _Diff _Count,
		_Wrap_alloc<_Alloc>& _Al, false_type)
	{	
	_FwdIt _Next = _First;

	try {
	for (; 0 < _Count; --_Count, (void)++_First)
		_Al.construct(_Unfancy(_First));
	} catch (...) {
	_Destroy_range(_Next, _First, _Al);
	throw;
	}
	}

template<class _FwdIt,
	class _Diff,
	class _Alloc> inline
	void _Uninitialized_default_fill_n1(_FwdIt _First, _Diff _Count,
		_Wrap_alloc<_Alloc>&, true_type)
	{	
	:: memset(_First, 0, _Count * sizeof(_Iter_value_t<_FwdIt>));
	}

template<class _FwdIt,
	class _Diff,
	class _Alloc> inline
	void _Uninitialized_default_fill_n(_FwdIt _First, _Diff _Count,
		_Wrap_alloc<_Alloc>& _Al)
	{	
	typedef _Iter_value_t<_FwdIt> _Ty;
	_Uninitialized_default_fill_n1(_First, _Count, _Al,
		typename conjunction<
			is_pointer<_FwdIt>,
			is_scalar<_Ty>,
			negation<is_volatile<_Ty>>,
			negation<is_member_pointer<_Ty>>,
			_Uses_default_construct<_Alloc, decltype(_Unfancy(_First))>>::type());
	}

		
template<class _OutIt,
	class _Ty>
	class raw_storage_iterator
		: public _Outit
	{	
public:
	explicit raw_storage_iterator(_OutIt _First)
		: _Next(_First)
		{	
		}

	raw_storage_iterator& operator*()
		{	
		return (*this);
		}

	raw_storage_iterator& operator=(const _Ty& _Val)
		{	
		_Construct(_Unfancy(_Next), _Val);
		return (*this);
		}

	raw_storage_iterator& operator=(_Ty&& _Val)
		{	
		_Construct(_Unfancy(_Next), ::std:: move(_Val));
		return (*this);
		}

	raw_storage_iterator& operator++()
		{	
		++_Next;
		return (*this);
		}

	raw_storage_iterator operator++(int)
		{	
		raw_storage_iterator _Ans = *this;
		++_Next;
		return (_Ans);
		}

	_OutIt base() const
		{	
		return (_Next);
		}

private:
	_OutIt _Next;	
	};

		
template<class _Ty>
	class _Temp_iterator
		: public _Outit
	{	
public:
	typedef _Ty *_Pty;

	_Temp_iterator(ptrdiff_t _Count = 0)
		{	
		_Buf._Begin = 0;
		_Buf._Current = 0;
		_Buf._Hiwater = 0;
		_Buf._Size = _Count;	
		_Pbuf = &_Buf;
		}

	_Temp_iterator(const _Temp_iterator& _Right)
		{	
		_Buf._Begin = 0;	
		_Buf._Current = 0;
		_Buf._Hiwater = 0;
		_Buf._Size = 0;
		*this = _Right;
		}

	~_Temp_iterator() noexcept
		{	
		if (_Buf._Begin != 0)
			{	
			for (_Pty _Next = _Buf._Begin;
				_Next != _Buf._Hiwater; ++_Next)
				_Destroy(_Next);
			::std:: return_temporary_buffer(_Buf._Begin);
			}
		}

	_Temp_iterator& operator=(const _Temp_iterator& _Right)
		{	
		_Pbuf = _Right._Pbuf;
		return (*this);
		}

	_Temp_iterator& operator=(const _Ty& _Val)
		{	
		if (_Pbuf->_Current < _Pbuf->_Hiwater)
			*_Pbuf->_Current++ = _Val;	
		else
			{	
			_Pty _Ptr = _Pbuf->_Current;
			_Construct(_Ptr, _Val);
			_Pbuf->_Hiwater = ++_Pbuf->_Current;
			}

		return (*this);
		}

	_Temp_iterator& operator=(_Ty&& _Val)
		{	
		if (_Pbuf->_Current < _Pbuf->_Hiwater)
			*_Pbuf->_Current++ =
				::std:: forward<_Ty>(_Val);	
		else
			{	
			_Pty _Ptr = _Pbuf->_Current;
			_Construct(_Ptr, ::std:: forward<_Ty>(_Val));
			_Pbuf->_Hiwater = ++_Pbuf->_Current;
			}

		return (*this);
		}

	_Temp_iterator& operator*()
		{	
		return (*this);
		}

	_Temp_iterator& operator++()
		{	
		return (*this);
		}

	_Temp_iterator& operator++(int)
		{	
		return (*this);
		}

	_Temp_iterator& _Init()
		{	
		_Pbuf->_Current = _Pbuf->_Begin;
		return (*this);
		}

	_Pty _First() const
		{	
		return (_Pbuf->_Begin);
		}

	_Pty _Last() const
		{	
		return (_Pbuf->_Current);
		}

	ptrdiff_t _Maxlen()
		{	
		if (_Pbuf->_Begin == 0 && 0 < _Pbuf->_Size)
			{	
			pair<_Pty, ptrdiff_t> _Pair =

				::std:: get_temporary_buffer<_Ty>(_Pbuf->_Size);

			_Pbuf->_Begin = _Pair.first;
			_Pbuf->_Current = _Pair.first;
			_Pbuf->_Hiwater = _Pair.first;
			_Pbuf->_Size = _Pair.second;
			}

		return (_Pbuf->_Size);
		}

private:
	struct _Bufpar
		{	
		_Pty _Begin;	
		_Pty _Current;	
		_Pty _Hiwater;	
		ptrdiff_t _Size;	
		};
	_Bufpar _Buf;	
	_Bufpar *_Pbuf;	
	};

 
		
template<class _Ty>
	class auto_ptr;

template<class _Ty>
	struct auto_ptr_ref
		{	
	explicit auto_ptr_ref(_Ty *_Right)
		: _Ref(_Right)
		{	
		}

	_Ty *_Ref;	
	};

template<class _Ty>
	class auto_ptr
		{	
public:
	typedef auto_ptr<_Ty> _Myt;
	typedef _Ty element_type;

	explicit auto_ptr(_Ty *_Ptr = 0) noexcept
		: _Myptr(_Ptr)
		{	
		}

	auto_ptr(_Myt& _Right) noexcept
		: _Myptr(_Right.release())
		{	
		}

	auto_ptr(auto_ptr_ref<_Ty> _Right) noexcept
		{	
		_Ty *_Ptr = _Right._Ref;
		_Right._Ref = 0;	
		_Myptr = _Ptr;	
		}

	template<class _Other>
		operator auto_ptr<_Other>() noexcept
		{	
		return (auto_ptr<_Other>(*this));
		}

	template<class _Other>
		operator auto_ptr_ref<_Other>() noexcept
		{	
		_Other *_Cvtptr = _Myptr;	
		auto_ptr_ref<_Other> _Ans(_Cvtptr);
		_Myptr = 0;	
		return (_Ans);
		}

	template<class _Other>
		_Myt& operator=(auto_ptr<_Other>& _Right) noexcept
		{	
		reset(_Right.release());
		return (*this);
		}

	template<class _Other>
		auto_ptr(auto_ptr<_Other>& _Right) noexcept
		: _Myptr(_Right.release())
		{	
		}

	_Myt& operator=(_Myt& _Right) noexcept
		{	
		reset(_Right.release());
		return (*this);
		}

	_Myt& operator=(auto_ptr_ref<_Ty> _Right) noexcept
		{	
		_Ty *_Ptr = _Right._Ref;
		_Right._Ref = 0;	
		reset(_Ptr);	
		return (*this);
		}

	~auto_ptr() noexcept
		{	
		delete _Myptr;
		}

	_Ty& operator*() const noexcept
		{	
 




		return (*get());
		}

	_Ty *operator->() const noexcept
		{	
 




		return (get());
		}

	_Ty *get() const noexcept
		{	
		return (_Myptr);
		}

	_Ty *release() noexcept
		{	
		_Ty *_Tmp = _Myptr;
		_Myptr = 0;
		return (_Tmp);
		}

	void reset(_Ty *_Ptr = 0)
		{	
		if (_Ptr != _Myptr)
			delete _Myptr;
		_Myptr = _Ptr;
		}

private:
	_Ty *_Myptr;	
	};
 
}

 
 #pragma warning(pop)
 #pragma pack(pop)










#pragma once






#pragma once





 #pragma pack(push,8)
 #pragma warning(push,3)
 
 

namespace std {
 #pragma warning(disable: 4127)

  #pragma warning(disable: 4251)

template<class _Elem,
	class _Traits = char_traits<_Elem>,
	class _Ax = allocator<_Elem> >
	class basic_string;

		
template<class _Mystr>
	class _String_const_iterator
		: public _Iterator012<random_access_iterator_tag,
			typename _Mystr::value_type,
			typename _Mystr::difference_type,
			typename _Mystr::const_pointer,
			typename _Mystr::const_reference,
			_Iterator_base>
	{	
public:
	typedef _String_const_iterator<_Mystr> _Myiter;
	typedef random_access_iterator_tag iterator_category;

	typedef typename _Mystr::value_type value_type;
	typedef typename _Mystr::difference_type difference_type;
	typedef typename _Mystr::const_pointer pointer;
	typedef typename _Mystr::const_reference reference;

	_String_const_iterator()
		: _Ptr()
		{	
		}

	_String_const_iterator(pointer _Parg, const _Container_base *_Pstring)
		: _Ptr(_Parg)
		{	
		this->_Adopt(_Pstring);
		}

	typedef pointer _Unchecked_type;

	_Myiter& _Rechecked(_Unchecked_type _Right)
		{	
		_Ptr = _Right;
		return (*this);
		}

	_Unchecked_type _Unchecked() const
		{	
		return (_Ptr);
		}

	reference operator*() const
		{	
 



















		;

		return (*_Ptr);
		}

	pointer operator->() const
		{	
		return (pointer_traits<pointer>::pointer_to(**this));
		}

	_Myiter& operator++()
		{	
 
















		++_Ptr;
		return (*this);
		}

	_Myiter operator++(int)
		{	
		_Myiter _Tmp = *this;
		++*this;
		return (_Tmp);
		}

	_Myiter& operator--()
		{	
 
















		--_Ptr;
		return (*this);
		}

	_Myiter operator--(int)
		{	
		_Myiter _Tmp = *this;
		--*this;
		return (_Tmp);
		}

	_Myiter& operator+=(difference_type _Off)
		{	
 


























		_Ptr += _Off;
		return (*this);
		}

	_Myiter operator+(difference_type _Off) const
		{	
		_Myiter _Tmp = *this;
		return (_Tmp += _Off);
		}

	_Myiter& operator-=(difference_type _Off)
		{	
		return (*this += -_Off);
		}

	_Myiter operator-(difference_type _Off) const
		{	
		_Myiter _Tmp = *this;
		return (_Tmp -= _Off);
		}

	difference_type operator-(const _Myiter& _Right) const
		{	
		_Compat(_Right);
		return (_Ptr - _Right._Ptr);
		}

	reference operator[](difference_type _Off) const
		{	
		return (*(*this + _Off));
		}

	bool operator==(const _Myiter& _Right) const
		{	
		_Compat(_Right);
		return (_Ptr == _Right._Ptr);
		}

	bool operator!=(const _Myiter& _Right) const
		{	
		return (!(*this == _Right));
		}

	bool operator<(const _Myiter& _Right) const
		{	
		_Compat(_Right);
		return (_Ptr < _Right._Ptr);
		}

	bool operator>(const _Myiter& _Right) const
		{	
		return (_Right < *this);
		}

	bool operator<=(const _Myiter& _Right) const
		{	
		return (!(_Right < *this));
		}

	bool operator>=(const _Myiter& _Right) const
		{	
		return (!(*this < _Right));
		}

 
















	void _Compat(const _Myiter&) const
		{	
		}
 

	pointer _Ptr;	
	};

template<class _Mystr> inline
	typename _String_const_iterator<_Mystr>::_Unchecked_type
		_Unchecked(_String_const_iterator<_Mystr> _Iter)
	{	
	return (_Iter._Unchecked());
	}

template<class _Mystr> inline
	_String_const_iterator<_Mystr>
		_Rechecked(_String_const_iterator<_Mystr>& _Iter,
			typename _String_const_iterator<_Mystr>
				::_Unchecked_type _Right)
	{	
	return (_Iter._Rechecked(_Right));
	}

template<class _Mystr> inline
	_String_const_iterator<_Mystr> operator+(
		typename _String_const_iterator<_Mystr>
			::difference_type _Off,
		_String_const_iterator<_Mystr> _Next)
	{	
	return (_Next += _Off);
	}

		
template<class _Mystr>
	class _String_iterator
		: public _String_const_iterator<_Mystr>
	{	
public:
	typedef _String_iterator<_Mystr> _Myiter;
	typedef _String_const_iterator<_Mystr> _Mybase;
	typedef random_access_iterator_tag iterator_category;

	typedef typename _Mystr::value_type value_type;
	typedef typename _Mystr::difference_type difference_type;
	typedef typename _Mystr::pointer pointer;
	typedef typename _Mystr::reference reference;

	_String_iterator()
		{	
		}

	_String_iterator(pointer _Parg, const _Container_base *_Pstring)
		: _Mybase(_Parg, _Pstring)
		{	
		}

	typedef pointer _Unchecked_type;

	_Myiter& _Rechecked(_Unchecked_type _Right)
		{	
		this->_Ptr = _Right;
		return (*this);
		}

	_Unchecked_type _Unchecked() const
		{	
		return (_Const_cast(this->_Ptr));
		}

	reference operator*() const
		{	
		return ((reference)**(_Mybase *)this);
		}

	pointer operator->() const
		{	
		return (pointer_traits<pointer>::pointer_to(**this));
		}

	_Myiter& operator++()
		{	
		++*(_Mybase *)this;
		return (*this);
		}

	_Myiter operator++(int)
		{	
		_Myiter _Tmp = *this;
		++*this;
		return (_Tmp);
		}

	_Myiter& operator--()
		{	
		--*(_Mybase *)this;
		return (*this);
		}

	_Myiter operator--(int)
		{	
		_Myiter _Tmp = *this;
		--*this;
		return (_Tmp);
		}

	_Myiter& operator+=(difference_type _Off)
		{	
		*(_Mybase *)this += _Off;
		return (*this);
		}

	_Myiter operator+(difference_type _Off) const
		{	
		_Myiter _Tmp = *this;
		return (_Tmp += _Off);
		}

	_Myiter& operator-=(difference_type _Off)
		{	
		return (*this += -_Off);
		}

	_Myiter operator-(difference_type _Off) const
		{	
		_Myiter _Tmp = *this;
		return (_Tmp -= _Off);
		}

	difference_type operator-(const _Mybase& _Right) const
		{	
		return (*(_Mybase *)this - _Right);
		}

	reference operator[](difference_type _Off) const
		{	
		return (*(*this + _Off));
		}
	};

template<class _Mystr> inline
	typename _String_iterator<_Mystr>::_Unchecked_type
		_Unchecked(_String_iterator<_Mystr> _Iter)
	{	
	return (_Iter._Unchecked());
	}

template<class _Mystr> inline
	_String_iterator<_Mystr>
		_Rechecked(_String_iterator<_Mystr>& _Iter,
			typename _String_iterator<_Mystr>
				::_Unchecked_type _Right)
	{	
	return (_Iter._Rechecked(_Right));
	}

template<class _Mystr> inline
	_String_iterator<_Mystr> operator+(
		typename _String_iterator<_Mystr>
			::difference_type _Off,
		_String_iterator<_Mystr> _Next)
	{	
	return (_Next += _Off);
	}

		
template<class _Value_type,
	class _Size_type,
	class _Difference_type,
	class _Pointer,
	class _Const_pointer,
	class _Reference,
	class _Const_reference>
	struct _String_iter_types
	{	
	typedef _Value_type value_type;
	typedef _Size_type size_type;
	typedef _Difference_type difference_type;
	typedef _Pointer pointer;
	typedef _Const_pointer const_pointer;
	typedef _Reference reference;
	typedef _Const_reference const_reference;
	};

template<class _Ty,
	class _Alloc0>
	struct _String_base_types
	{	
	typedef _Alloc0 _Alloc;
	typedef _String_base_types<_Ty, _Alloc> _Myt;

	typedef _Wrap_alloc<_Alloc> _Alty0;
	typedef typename _Alty0::template rebind<_Ty>::other _Alty;


	typedef typename _If<_Is_simple_alloc<_Alty>::value,
		_Simple_types<typename _Alty::value_type>,
		_String_iter_types<typename _Alty::value_type,
			typename _Alty::size_type,
			typename _Alty::difference_type,
			typename _Alty::pointer,
			typename _Alty::const_pointer,
			typename _Alty::reference,
			typename _Alty::const_reference> >::type
		_Val_types;
	};

		
template<class _Val_types>
	class _String_val
		: public _Container_base
	{	
public:
	typedef _String_val<_Val_types> _Myt;

	typedef typename _Val_types::value_type value_type;
	typedef typename _Val_types::size_type size_type;
	typedef typename _Val_types::difference_type difference_type;
	typedef typename _Val_types::pointer pointer;
	typedef typename _Val_types::const_pointer const_pointer;
	typedef typename _Val_types::reference reference;
	typedef typename _Val_types::const_reference const_reference;

	typedef _String_iterator<_Myt> iterator;
	typedef _String_const_iterator<_Myt> const_iterator;

	_String_val()
		: _Bx(),
		_Mysize(0),
		_Myres(0)
		{	
		}

	enum
		{	
		_BUF_SIZE = 16 / sizeof (value_type) < 1 ? 1
			: 16 / sizeof (value_type)};
	enum
		{	
		_ALLOC_MASK = sizeof (value_type) <= 1 ? 15
			: sizeof (value_type) <= 2 ? 7
			: sizeof (value_type) <= 4 ? 3
			: sizeof (value_type) <= 8 ? 1 : 0};

	value_type *_Myptr()
		{	
		return (this->_BUF_SIZE <= _Myres
			? _Unfancy(_Bx._Ptr)
			: _Bx._Buf);
		}

	const value_type *_Myptr() const
		{	
		return (this->_BUF_SIZE <= _Myres
			? _Unfancy(_Bx._Ptr)
			: _Bx._Buf);
		}

	union _Bxty
		{	
		_Bxty()
			{	
			}

		~_Bxty() noexcept
			{	
			}

		value_type _Buf[_BUF_SIZE];
		pointer _Ptr;
		char _Alias[_BUF_SIZE];	
		} _Bx;

	size_type _Mysize;	
	size_type _Myres;	
	};

		
template<class _Alloc_types>
	class _String_alloc
	{	
public:
	typedef _String_alloc<_Alloc_types> _Myt;
	typedef typename _Alloc_types::_Alloc _Alloc;
	typedef typename _Alloc_types::_Alty _Alty;
	typedef typename _Alloc_types::_Val_types _Val_types;

	typedef typename _Val_types::value_type value_type;
	typedef typename _Val_types::size_type size_type;
	typedef typename _Val_types::difference_type difference_type;
	typedef typename _Val_types::pointer pointer;
	typedef typename _Val_types::const_pointer const_pointer;
	typedef typename _Val_types::reference reference;
	typedef typename _Val_types::const_reference const_reference;

	typedef _String_iterator<_String_val<_Val_types> > iterator;
	typedef _String_const_iterator<_String_val<_Val_types> > const_iterator;

	enum
		{	
		_BUF_SIZE = _String_val<_Val_types>::_BUF_SIZE
		};

	enum
		{	
		_ALLOC_MASK = _String_val<_Val_types>::_ALLOC_MASK
		};

	value_type *_Myptr()
		{	
		return (_Get_data()._Myptr());
		}

	const value_type *_Myptr() const
		{	
		return (_Get_data()._Myptr());
		}

 
	_String_alloc()
		: _Mypair(_Zero_then_variadic_args_t())
		{	
		}

	template<class _Any_alloc,
		class = enable_if_t<!is_same<decay_t<_Any_alloc>, _Myt>::value> >
		_String_alloc(_Any_alloc&& _Al)
		: _Mypair(_One_then_variadic_args_t(),
			::std:: forward<_Any_alloc>(_Al))
		{	
		}

	void _Copy_alloc(const _Alty& _Al)
		{	
		_Pocca(_Getal(), _Al);
		}

	void _Move_alloc(_Alty& _Al)
		{	
		_Pocma(_Getal(), _Al);
		}

 





































































	void _Orphan_all()
		{	
		_Get_data()._Orphan_all();
		}

	void _Swap_all(_Myt& _Right)
		{	
		_Get_data()._Swap_all(_Right._Get_data());
		}

	_Alty& _Getal() noexcept
		{	
		return (_Mypair._Get_first());
		}

	const _Alty& _Getal() const noexcept
		{	
		return (_Mypair._Get_first());
		}

	_String_val<_Val_types>& _Get_data() noexcept
		{	
		return (_Mypair._Get_second());
		}

	const _String_val<_Val_types>& _Get_data() const noexcept
		{	
		return (_Mypair._Get_second());
		}

	typedef typename _String_val<_Val_types>::_Bxty _Bxty;

	_Bxty& _Bx() noexcept
		{	
		return (_Get_data()._Bx);
		}

	const _Bxty& _Bx() const noexcept
		{	
		return (_Get_data()._Bx);
		}

	size_type& _Mysize() noexcept
		{	
		return (_Get_data()._Mysize);
		}

	const size_type& _Mysize() const noexcept
		{	
		return (_Get_data()._Mysize);
		}

	size_type& _Myres() noexcept
		{	
		return (_Get_data()._Myres);
		}

	const size_type& _Myres() const noexcept
		{	
		return (_Get_data()._Myres);
		}

private:
	_Compressed_pair<_Alty, _String_val<_Val_types> > _Mypair;
	};

		
template<class _Elem,
	class _Traits,
	class _Alloc>
	class basic_string
		: public _String_alloc<_String_base_types<_Elem, _Alloc> >
	{	
public:
	typedef basic_string<_Elem, _Traits, _Alloc> _Myt;
	typedef _String_alloc<_String_base_types<_Elem, _Alloc> > _Mybase;
	typedef _Traits traits_type;
	typedef _Alloc allocator_type;

	typedef typename _Mybase::_Alty _Alty;

	typedef typename _Mybase::value_type value_type;
	typedef typename _Mybase::size_type size_type;
	typedef typename _Mybase::difference_type difference_type;
	typedef typename _Mybase::pointer pointer;
	typedef typename _Mybase::const_pointer const_pointer;
	typedef typename _Mybase::reference reference;
	typedef typename _Mybase::const_reference const_reference;

	typedef typename _Mybase::iterator iterator;
	typedef typename _Mybase::const_iterator const_iterator;

	typedef ::std:: reverse_iterator<iterator> reverse_iterator;
	typedef ::std:: reverse_iterator<const_iterator> const_reverse_iterator;

	basic_string(const _Myt& _Right)

		: _Mybase(_Right._Getal().select_on_container_copy_construction())


		{	
		_Tidy();
		assign(_Right, 0, npos);
		}

	basic_string(const _Myt& _Right, const _Alloc& _Al)
		: _Mybase(_Al)
		{	
		_Tidy();
		assign(_Right, 0, npos);
		}

	basic_string() noexcept(is_nothrow_default_constructible<_Alloc>::value)
		: _Mybase()
		{	
		_Tidy();
		}

	explicit basic_string(const _Alloc& _Al) noexcept
		: _Mybase(_Al)
		{	
		_Tidy();
		}

	basic_string(const _Myt& _Right, size_type _Roff,
		size_type _Count = npos)
		: _Mybase(_Right._Getal())
		{	
		_Tidy();
		assign(_Right, _Roff, _Count);
		}

	basic_string(const _Myt& _Right, size_type _Roff, size_type _Count,
		const _Alloc& _Al)
		: _Mybase(_Al)
		{	
		_Tidy();
		assign(_Right, _Roff, _Count);
		}

	basic_string(const _Elem *_Ptr, size_type _Count)
		: _Mybase()
		{	
		_Tidy();
		assign(_Ptr, _Count);
		}

	basic_string(const _Elem *_Ptr, size_type _Count, const _Alloc& _Al)
		: _Mybase(_Al)
		{	
		_Tidy();
		assign(_Ptr, _Count);
		}

	basic_string(const _Elem *_Ptr)
		: _Mybase()
		{	
		_Tidy();
		assign(_Ptr);
		}

	basic_string(const _Elem *_Ptr, const _Alloc& _Al)
		: _Mybase(_Al)
		{	
		_Tidy();
		assign(_Ptr);
		}

	basic_string(size_type _Count, _Elem _Ch)
		: _Mybase()
		{	
		_Tidy();
		assign(_Count, _Ch);
		}

	basic_string(size_type _Count, _Elem _Ch, const _Alloc& _Al)
		: _Mybase(_Al)
		{	
		_Tidy();
		assign(_Count, _Ch);
		}

	template<class _Iter,
		class = typename enable_if<_Is_iterator<_Iter>::value,
			void>::type>
		basic_string(_Iter _First, _Iter _Last, const _Alloc& _Al = _Alloc())
		: _Mybase(_Al)
		{	
		;
		_Tidy();
		_Construct(_Unchecked(_First), _Unchecked(_Last), _Iter_cat_t<_Iter>());
		}

	template<class _Iter>
		void _Construct(_Iter _First,
			_Iter _Last, input_iterator_tag)
		{	
		try {
		for (; _First != _Last; ++_First)
			append((size_type)1, (_Elem)*_First);
		} catch (...) {
		_Tidy(true);
		throw;
		}
		}

	template<class _Iter>
		void _Construct(_Iter _First,
			_Iter _Last, forward_iterator_tag)
		{	
		size_type _Count = ::std:: distance(_First, _Last);
		reserve(_Count);
		_Construct(_First, _Last, input_iterator_tag());
		}

	void _Construct(_Elem *_First,
		_Elem *_Last, random_access_iterator_tag)
		{	
		if (_First != _Last)
			assign(_First, _Last - _First);
		}

	void _Construct(const _Elem *_First,
		const _Elem *_Last, random_access_iterator_tag)
		{	
		if (_First != _Last)
			assign(_First, _Last - _First);
		}

	basic_string(_Myt&& _Right) noexcept
		: _Mybase(::std:: move(_Right._Getal()))
		{	
		_Tidy();
		_Assign_rv(::std:: forward<_Myt>(_Right));
		}

	basic_string(_Myt&& _Right, const _Alloc& _Al)
		: _Mybase(_Al)
		{	
		if (this->_Getal() != _Right._Getal())
			assign(_Right.begin(), _Right.end());
		else
			_Assign_rv(::std:: forward<_Myt>(_Right));
		}

	_Myt& operator=(_Myt&& _Right)
		noexcept(_Alty::propagate_on_container_move_assignment::value || _Alty::is_always_equal::value)
		{	
		if (this != &_Right)
			{	
			_Tidy(true);

			if (_Alty::propagate_on_container_move_assignment::value
				&& this->_Getal() != _Right._Getal())
				this->_Move_alloc(_Right._Getal());

			if (this->_Getal() != _Right._Getal())
				assign(_Right.begin(), _Right.end());
			else
				_Assign_rv(::std:: forward<_Myt>(_Right));
			}
		return (*this);
		}

	_Myt& assign(_Myt&& _Right) noexcept
		{	
		if (this == &_Right)
			;
		else if (get_allocator() != _Right.get_allocator()
			&& this->_BUF_SIZE <= _Right._Myres())
			*this = _Right;
		else
			{	
			_Tidy(true);
			_Assign_rv(::std:: forward<_Myt>(_Right));
			}
		return (*this);
		}

	void _Assign_rv(_Myt&& _Right)
		{	
		if (_Right._Myres() < this->_BUF_SIZE)
			_Traits::move(this->_Bx()._Buf, _Right._Bx()._Buf,
				_Right._Mysize() + 1);
		else
			{	
			this->_Getal().construct(::std:: addressof(this->_Bx()._Ptr), _Right._Bx()._Ptr);
			_Right._Bx()._Ptr = pointer();
			}
		this->_Mysize() = _Right._Mysize();
		this->_Myres() = _Right._Myres();
		_Right._Tidy();
		}

	basic_string(::std:: initializer_list<_Elem> _Ilist,
		const _Alloc& _Al = allocator_type())
		: _Mybase(_Al)
		{	
		_Tidy();
		assign(_Ilist.begin(), _Ilist.end());
		}

	_Myt& operator=(::std:: initializer_list<_Elem> _Ilist)
		{	
		return (assign(_Ilist.begin(), _Ilist.end()));
		}

	_Myt& operator+=(::std:: initializer_list<_Elem> _Ilist)
		{	
		return (append(_Ilist.begin(), _Ilist.end()));
		}

	_Myt& assign(::std:: initializer_list<_Elem> _Ilist)
		{	
		return (assign(_Ilist.begin(), _Ilist.end()));
		}

	_Myt& append(::std:: initializer_list<_Elem> _Ilist)
		{	
		return (append(_Ilist.begin(), _Ilist.end()));
		}

	iterator insert(const_iterator _Where,
		::std:: initializer_list<_Elem> _Ilist)
		{	
		return (insert(_Where, _Ilist.begin(), _Ilist.end()));
		}

	_Myt& replace(const_iterator _First, const_iterator _Last,
		::std:: initializer_list<_Elem> _Ilist)
		{	
		return (replace(_First, _Last, _Ilist.begin(), _Ilist.end()));
		}

	~basic_string() noexcept
		{	
		_Tidy(true);
		}

	 static const size_type npos;	

	_Myt& operator=(const _Myt& _Right)
		{	
		if (this != &_Right)
			{	
			if (this->_Getal() != _Right._Getal()
				&& _Alty::propagate_on_container_copy_assignment::value)
				{	
				_Tidy(true);
				this->_Copy_alloc(_Right._Getal());
				}

			assign(_Right);
			}
		return (*this);
		}

	_Myt& operator=(const _Elem *_Ptr)
		{	
		return (assign(_Ptr));
		}

	_Myt& operator=(_Elem _Ch)
		{	
		return (assign(1, _Ch));
		}

	_Myt& operator+=(const _Myt& _Right)
		{	
		return (append(_Right));
		}

	_Myt& operator+=(const _Elem *_Ptr)
		{	
		return (append(_Ptr));
		}

	_Myt& operator+=(_Elem _Ch)
		{	
		return (append((size_type)1, _Ch));
		}

	_Myt& append(const _Myt& _Right)
		{	
		return (append(_Right, 0, npos));
		}

	_Myt& append(const _Myt& _Right,
		size_type _Roff, size_type _Count = npos)
		{	
		_Right._Check_offset(_Roff);
		_Count = _Right._Clamp_suffix_size(_Roff, _Count);
		if (npos - this->_Mysize() <= _Count)
			_Xlen();	

		const size_type _Num = this->_Mysize() + _Count;
		if (0 < _Count && _Grow(_Num))
			{	
			_Traits::copy(this->_Myptr() + this->_Mysize(),
				_Right._Myptr() + _Roff, _Count);
			_Eos(_Num);
			}
		return (*this);
		}

	_Myt& append(const _Elem *_Ptr, size_type _Count)
		{	
		;
		if (_Inside(_Ptr))
			return (append(*this,
				_Ptr - this->_Myptr(), _Count));	
		if (npos - this->_Mysize() <= _Count)
			_Xlen();	

		const size_type _Num = this->_Mysize() + _Count;
		if (0 < _Count && _Grow(_Num))
			{	
			_Traits::copy(this->_Myptr() + this->_Mysize(), _Ptr, _Count);
			_Eos(_Num);
			}
		return (*this);
		}

	_Myt& append(const _Elem *_Ptr)
		{	
		;
		return (append(_Ptr, _Traits::length(_Ptr)));
		}

	_Myt& append(size_type _Count, _Elem _Ch)
		{	
		if (npos - this->_Mysize() <= _Count)
			_Xlen();	

		const size_type _Num = this->_Mysize() + _Count;
		if (0 < _Count && _Grow(_Num))
			{	
			_Chassign(this->_Mysize(), _Count, _Ch);
			_Eos(_Num);
			}
		return (*this);
		}

	template<class _Iter>
		typename enable_if<_Is_iterator<_Iter>::value,
			_Myt&>::type
		append(_Iter _First, _Iter _Last)
		{	
		return (replace(end(), end(), _First, _Last));
		}

	_Myt& append(const_pointer _First, const_pointer _Last)
		{	
		return (replace(end(), end(), _First, _Last));
		}

	_Myt& append(const_iterator _First, const_iterator _Last)
		{	
		return (replace(end(), end(), _First, _Last));
		}

	_Myt& assign(const _Myt& _Right)
		{	
		return (assign(_Right, 0, npos));
		}

	_Myt& assign(const _Myt& _Right,
		size_type _Roff, size_type _Count = npos)
		{	
		_Right._Check_offset(_Roff);
		_Count = _Right._Clamp_suffix_size(_Roff, _Count);

		if (this == &_Right)
			erase((size_type)(_Roff + _Count)), erase(0, _Roff);	
		else if (_Grow(_Count))
			{	
			_Traits::copy(this->_Myptr(),
				_Right._Myptr() + _Roff, _Count);
			_Eos(_Count);
			}
		return (*this);
		}

	_Myt& assign(const _Elem *_Ptr, size_type _Count)
		{	
		;
		if (_Inside(_Ptr))
			return (assign(*this,
				_Ptr - this->_Myptr(), _Count));	

		if (_Grow(_Count))
			{	
			_Traits::copy(this->_Myptr(), _Ptr, _Count);
			_Eos(_Count);
			}
		return (*this);
		}

	_Myt& assign(const _Elem *_Ptr)
		{	
		;
		return (assign(_Ptr, _Traits::length(_Ptr)));
		}

	_Myt& assign(size_type _Count, _Elem _Ch)
		{	
		if (_Count == npos)
			_Xlen();	

		if (_Grow(_Count))
			{	
			_Chassign(0, _Count, _Ch);
			_Eos(_Count);
			}
		return (*this);
		}

	template<class _Iter>
		typename enable_if<_Is_iterator<_Iter>::value,
			_Myt&>::type
		assign(_Iter _First, _Iter _Last)
		{	
		return (replace(begin(), end(), _First, _Last));
		}

	_Myt& assign(const_pointer _First, const_pointer _Last)
		{	
		return (replace(begin(), end(), _First, _Last));
		}

	_Myt& assign(const_iterator _First, const_iterator _Last)
		{	
		return (replace(begin(), end(), _First, _Last));
		}

	_Myt& insert(size_type _Off, const _Myt& _Right)
		{	
		return (insert(_Off, _Right, 0, npos));
		}

	_Myt& insert(size_type _Off,
		const _Myt& _Right, size_type _Roff, size_type _Count = npos)
		{	
		_Check_offset(_Off);
		_Right._Check_offset(_Roff);
		_Count = _Right._Clamp_suffix_size(_Roff, _Count);
		if (npos - this->_Mysize() <= _Count)
			_Xlen();	

		const size_type _Num = this->_Mysize() + _Count;
		if (0 < _Count && _Grow(_Num))
			{	
			_Traits::move(this->_Myptr() + _Off + _Count,
				this->_Myptr() + _Off,
				this->_Mysize() - _Off);	
			if (this == &_Right)
				_Traits::move(this->_Myptr() + _Off,
					this->_Myptr() + (_Off < _Roff ? _Roff + _Count : _Roff),
						_Count);	
			else
				_Traits::copy(this->_Myptr() + _Off,
					_Right._Myptr() + _Roff, _Count);	
			_Eos(_Num);
			}
		return (*this);
		}

	_Myt& insert(size_type _Off,
		const _Elem *_Ptr, size_type _Count)
		{	
		;
		if (_Inside(_Ptr))
			return (insert(_Off, *this,
				_Ptr - this->_Myptr(), _Count));	
		_Check_offset(_Off);
		if (npos - this->_Mysize() <= _Count)
			_Xlen();	
		const size_type _Num = this->_Mysize() + _Count;
		if (0 < _Count && _Grow(_Num))
			{	
			_Traits::move(this->_Myptr() + _Off + _Count,
				this->_Myptr() + _Off,
				this->_Mysize() - _Off);	
			_Traits::copy(this->_Myptr() + _Off, _Ptr, _Count);	
			_Eos(_Num);
			}
		return (*this);
		}

	_Myt& insert(size_type _Off, const _Elem *_Ptr)
		{	
		;
		return (insert(_Off, _Ptr, _Traits::length(_Ptr)));
		}

	_Myt& insert(size_type _Off,
		size_type _Count, _Elem _Ch)
		{	
		_Check_offset(_Off);
		if (npos - this->_Mysize() <= _Count)
			_Xlen();	
		const size_type _Num = this->_Mysize() + _Count;
		if (0 < _Count && _Grow(_Num))
			{	
			_Traits::move(this->_Myptr() + _Off + _Count,
				this->_Myptr() + _Off,
				this->_Mysize() - _Off);	
			_Chassign(_Off, _Count, _Ch);	
			_Eos(_Num);
			}
		return (*this);
		}

	iterator insert(const_iterator _Where)
		{	
		return (insert(_Where, _Elem()));
		}

	iterator insert(const_iterator _Where, _Elem _Ch)
		{	
		size_type _Off = _Where - begin();
		insert(_Off, 1, _Ch);
		return (begin() + _Off);
		}

	iterator insert(const_iterator _Where, size_type _Count, _Elem _Ch)
		{	
		size_type _Off = _Where - begin();
		insert(_Off, _Count, _Ch);
		return (begin() + _Off);
		}

	template<class _Iter>
		typename enable_if<_Is_iterator<_Iter>::value,
			iterator>::type
		insert(const_iterator _Where, _Iter _First, _Iter _Last)
		{	
		size_type _Off = _Where - begin();
		replace(_Where, _Where, _First, _Last);
		return (begin() + _Off);
		}

	iterator insert(const_iterator _Where,
		const_pointer _First, const_pointer _Last)
		{	
		size_type _Off = _Where - begin();
		replace(_Where, _Where, _First, _Last);
		return (begin() + _Off);
		}

	iterator insert(const_iterator _Where,
		const_iterator _First, const_iterator _Last)
		{	
		size_type _Off = _Where - begin();
		replace(_Where, _Where, _First, _Last);
		return (begin() + _Off);
		}

	_Myt& erase(size_type _Off = 0)
		{	
		_Check_offset(_Off);
		_Eos(_Off);
		return (*this);
		}

	_Myt& erase(size_type _Off, size_type _Count)
		{	
		_Check_offset(_Off);
		if (this->_Mysize() - _Off <= _Count)
			_Eos(_Off);	
		else if (0 < _Count)
			{	
			value_type *_Ptr = this->_Myptr() + _Off;
			size_type _Newsize = this->_Mysize() - _Count;
			_Traits::move(_Ptr, _Ptr + _Count, _Newsize - _Off);
			_Eos(_Newsize);
			}
		return (*this);
		}

	iterator erase(const_iterator _Where)
		{	
		size_type _Count = _Where - begin();
		erase(_Count, 1);
		return (begin() + _Count);
		}

	iterator erase(const_iterator _First, const_iterator _Last)
		{	
		;
		size_type _Count = _First - begin();
		erase(_Count, _Last - _First);
		return (begin() + _Count);
		}

	void clear() noexcept
		{	
		_Eos(0);
		}

	_Myt& replace(size_type _Off, size_type _N0, const _Myt& _Right)
		{	
		return (replace(_Off, _N0, _Right, 0, npos));
		}

	_Myt& replace(size_type _Off,
		size_type _N0, const _Myt& _Right, size_type _Roff,
			size_type _Count = npos)
		{	
		_Check_offset(_Off);
		_Right._Check_offset(_Roff);
		_N0 = _Clamp_suffix_size(_Off, _N0);
		_Count = _Right._Clamp_suffix_size(_Roff, _Count);
		if (npos - _Count <= this->_Mysize() - _N0)
			_Xlen();	

		const size_type _Nm = this->_Mysize() - _N0 - _Off;	
		const size_type _Newsize = this->_Mysize() + _Count - _N0;
		if (this->_Mysize() < _Newsize)
			_Grow(_Newsize);

		if (_Count == _N0)
			{	
			_Traits::move(this->_Myptr() + _Off,
				_Right._Myptr() + _Roff, _Count);	
			}
		else if (this != &_Right)
			{	
			_Traits::move(this->_Myptr() + _Off + _Count,
				this->_Myptr() + _Off + _N0, _Nm);	
			_Traits::copy(this->_Myptr() + _Off,
				_Right._Myptr() + _Roff, _Count);	
			}
		else if (_Count < _N0)
			{	
			_Traits::move(this->_Myptr() + _Off,
				this->_Myptr() + _Roff, _Count);	
			_Traits::move(this->_Myptr() + _Off + _Count,
				this->_Myptr() + _Off + _N0, _Nm);	
			}
		else if (_Roff <= _Off)
			{	
			_Traits::move(this->_Myptr() + _Off + _Count,
				this->_Myptr() + _Off + _N0, _Nm);	
			_Traits::move(this->_Myptr() + _Off,
				this->_Myptr() + _Roff, _Count);	
			}
		else if (_Off + _N0 <= _Roff)
			{	
			_Traits::move(this->_Myptr() + _Off + _Count,
				this->_Myptr() + _Off + _N0, _Nm);	
			_Traits::move(this->_Myptr() + _Off,
				this->_Myptr() + (_Roff + _Count - _N0),
				_Count);	
			}
		else
			{	
			_Traits::move(this->_Myptr() + _Off,
				this->_Myptr() + _Roff, _N0);	
			_Traits::move(this->_Myptr() + _Off + _Count,
				this->_Myptr() + _Off + _N0, _Nm);	
			_Traits::move(this->_Myptr() + _Off + _N0,
				this->_Myptr() + _Roff + _Count,
				_Count - _N0);	
			}

		_Eos(_Newsize);
		return (*this);
		}

	_Myt& replace(size_type _Off,
		size_type _N0, const _Elem *_Ptr, size_type _Count)
		{	
		;
		if (_Inside(_Ptr))
			return (replace(_Off, _N0, *this,
				_Ptr - this->_Myptr(),
				_Count));	
		_Check_offset(_Off);
		_N0 = _Clamp_suffix_size(_Off, _N0);
		if (npos - _Count <= this->_Mysize() - _N0)
			_Xlen();	
		size_type _Nm = this->_Mysize() - _N0 - _Off;

		if (_Count < _N0)
			_Traits::move(this->_Myptr() + _Off + _Count,
				this->_Myptr() + _Off + _N0,
				_Nm);	
		const size_type _Num = this->_Mysize() + _Count - _N0;
		if ((0 < _Count || 0 < _N0)
			&& _Grow(_Num))
			{	
			if (_N0 < _Count)
				_Traits::move(this->_Myptr() + _Off + _Count,
					this->_Myptr() + _Off + _N0, _Nm);	
			_Traits::copy(this->_Myptr() + _Off, _Ptr, _Count);	
			_Eos(_Num);
			}
		return (*this);
		}

	_Myt& replace(size_type _Off, size_type _N0, const _Elem *_Ptr)
		{	
		;
		return (replace(_Off, _N0, _Ptr, _Traits::length(_Ptr)));
		}

	_Myt& replace(size_type _Off,
		size_type _N0, size_type _Count, _Elem _Ch)
		{	
		_Check_offset(_Off);
		_N0 = _Clamp_suffix_size(_Off, _N0);
		if (npos - _Count <= this->_Mysize() - _N0)
			_Xlen();	
		size_type _Nm = this->_Mysize() - _N0 - _Off;

		if (_Count < _N0)
			_Traits::move(this->_Myptr() + _Off + _Count,
				this->_Myptr() + _Off + _N0,
				_Nm);	
		const size_type _Num = this->_Mysize() + _Count - _N0;
		if ((0 < _Count || 0 < _N0)
			&& _Grow(_Num))
			{	
			if (_N0 < _Count)
				_Traits::move(this->_Myptr() + _Off + _Count,
					this->_Myptr() + _Off + _N0, _Nm);	
			_Chassign(_Off, _Count, _Ch);	
			_Eos(_Num);
			}
		return (*this);
		}

	_Myt& replace(const_iterator _First, const_iterator _Last,
		const _Myt& _Right)
		{	
		return (replace(_First - begin(), _Last - _First, _Right));
		}

	_Myt& replace(const_iterator _First, const_iterator _Last,
		const _Elem *_Ptr, size_type _Count)
		{	
		return (replace(_First - begin(), _Last - _First, _Ptr, _Count));
		}

	_Myt& replace(const_iterator _First, const_iterator _Last,
		const _Elem *_Ptr)
		{	
		return (replace(_First - begin(), _Last - _First, _Ptr));
		}

	_Myt& replace(const_iterator _First, const_iterator _Last,
		size_type _Count, _Elem _Ch)
		{	
		return (replace(_First - begin(), _Last - _First, _Count, _Ch));
		}

	template<class _Iter>
		typename enable_if<_Is_iterator<_Iter>::value,
			_Myt&>::type
		replace(const_iterator _First, const_iterator _Last,
			_Iter _First2, _Iter _Last2)
		{	
		_Myt _Right(_First2, _Last2);
		replace(_First, _Last, _Right);
		return (*this);
		}

	_Myt& replace(const_iterator _First, const_iterator _Last,
		const_pointer _First2, const_pointer _Last2)
		{	
		if (_First2 == _Last2)
			erase(_First - begin(), _Last - _First);
		else
			replace(_First - begin(), _Last - _First,
				&*_First2, _Last2 - _First2);
		return (*this);
		}

	_Myt& replace(const_iterator _First, const_iterator _Last,
		pointer _First2, pointer _Last2)
		{	
		if (_First2 == _Last2)
			erase(_First - begin(), _Last - _First);
		else
			replace(_First - begin(), _Last - _First,
				&*_First2, _Last2 - _First2);
		return (*this);
		}

	_Myt& replace(const_iterator _First, const_iterator _Last,
		const_iterator _First2, const_iterator _Last2)
		{	
		if (_First2 == _Last2)
			erase(_First - begin(), _Last - _First);
		else
			replace(_First - begin(), _Last - _First,
				&*_First2, _Last2 - _First2);
		return (*this);
		}

	_Myt& replace(const_iterator _First, const_iterator _Last,
		iterator _First2, iterator _Last2)
		{	
		if (_First2 == _Last2)
			erase(_First - begin(), _Last - _First);
		else
			replace(_First - begin(), _Last - _First,
				&*_First2, _Last2 - _First2);
		return (*this);
		}

	iterator begin() noexcept
		{	
		auto _Mydata = &this->_Get_data();
		return (iterator(this->_Getal().address(*_Mydata->_Myptr()), _Mydata));
		}

	const_iterator begin() const noexcept
		{	
		auto _Mydata = &this->_Get_data();
		return (const_iterator(this->_Getal().address(*_Mydata->_Myptr()), _Mydata));
		}

	iterator end() noexcept
		{	
		auto _Mydata = &this->_Get_data();
		return (iterator(this->_Getal().address(*_Mydata->_Myptr()) + _Mydata->_Mysize, _Mydata));
		}

	const_iterator end() const noexcept
		{	
		auto _Mydata = &this->_Get_data();
		return (const_iterator(this->_Getal().address(*_Mydata->_Myptr()) + _Mydata->_Mysize, _Mydata));
		}

	reverse_iterator rbegin() noexcept
		{	
		return (reverse_iterator(end()));
		}

	const_reverse_iterator rbegin() const noexcept
		{	
		return (const_reverse_iterator(end()));
		}

	reverse_iterator rend() noexcept
		{	
		return (reverse_iterator(begin()));
		}

	const_reverse_iterator rend() const noexcept
		{	
		return (const_reverse_iterator(begin()));
		}

	const_iterator cbegin() const noexcept
		{	
		return (begin());
		}

	const_iterator cend() const noexcept
		{	
		return (end());
		}

	const_reverse_iterator crbegin() const noexcept
		{	
		return (rbegin());
		}

	const_reverse_iterator crend() const noexcept
		{	
		return (rend());
		}

	void shrink_to_fit()
		{	
		if ((size() | this->_ALLOC_MASK) < capacity())
			{	
			_Myt _Tmp(*this);
			swap(_Tmp);
			}
		}

	reference at(size_type _Off)
		{	
		_Check_offset_exclusive(_Off);
		return (this->_Myptr()[_Off]);
		}

	const_reference at(size_type _Off) const
		{	
		_Check_offset_exclusive(_Off);
		return (this->_Myptr()[_Off]);
		}

	reference operator[](size_type _Off)
		{	
 







		return (this->_Myptr()[_Off]);
		}

	const_reference operator[](size_type _Off) const
		{	
 







		return (this->_Myptr()[_Off]);
		}

	void push_back(_Elem _Ch)
		{	
		auto& _Dx = this->_Get_data();
		auto& _Sz = _Dx._Mysize;
		if (_Sz == _Dx._Myres)
			_Grow(_Sz + 1); 
		auto _Ptr = _Dx._Myptr();
		_Traits::assign(_Ptr[_Sz], _Ch);
		++_Sz;
		_Traits::assign(_Ptr[_Sz], _Elem());
		}

	void pop_back()
		{	
		erase(this->_Mysize() - 1);	
		}

	reference front()
		{	
		return (*begin());
		}

	const_reference front() const
		{	
		return (*begin());
		}

	reference back()
		{	
		return (*(end() - 1));
		}

	const_reference back() const
		{	
		return (*(end() - 1));
		}

	const _Elem *c_str() const noexcept
		{	
		return (this->_Myptr());
		}

	const _Elem *data() const noexcept
		{	
		return (this->_Myptr());
		}








	size_type length() const noexcept
		{	
		return (this->_Mysize());
		}

	size_type size() const noexcept
		{	
		return (this->_Mysize());
		}

	size_type max_size() const noexcept
		{	
		const size_type _Num = this->_Getal().max_size();
		return (_Num <= 1 ? 1 : _Num - 1);
		}

	void resize(size_type _Newsize)
		{	
		resize(_Newsize, _Elem());
		}

	void resize(size_type _Newsize, _Elem _Ch)
		{	
		if (_Newsize <= this->_Mysize())
			_Eos(_Newsize);
		else
			append(_Newsize - this->_Mysize(), _Ch);
		}

	size_type capacity() const noexcept
		{	
		return (this->_Myres());
		}

	void reserve(size_type _Newcap = 0)
		{	
		if (this->_Mysize() <= _Newcap && this->_Myres() != _Newcap)
			{	
			size_type _Size = this->_Mysize();
			if (_Grow(_Newcap, true))
				_Eos(_Size);
			}
		}

	bool empty() const noexcept
		{	
		return (this->_Mysize() == 0);
		}

	
	size_type copy(_Elem *_Ptr,
		size_type _Count, size_type _Off = 0) const
		{	
		;
		_Check_offset(_Off);
		_Count = _Clamp_suffix_size(_Off, _Count);
		_Traits::copy(_Ptr, this->_Myptr() + _Off, _Count);
		return (_Count);
		}

	size_type _Copy_s(_Elem *_Dest, size_type _Dest_size,
		size_type _Count, size_type _Off = 0) const
		{	
		;
		_Check_offset(_Off);
		_Count = _Clamp_suffix_size(_Off, _Count);
		_Traits::_Copy_s(_Dest, _Dest_size, this->_Myptr() + _Off, _Count);
		return (_Count);
		}

	void _Swap_bx(_Myt& _Right)
		{	
		if (this->_BUF_SIZE <= this->_Myres())
			if (this->_BUF_SIZE <= _Right._Myres())
				_Swap_adl(this->_Bx()._Ptr, _Right._Bx()._Ptr);
			else
				{	
				pointer _Ptr = this->_Bx()._Ptr;
				this->_Getal().destroy(::std:: addressof(this->_Bx()._Ptr));
				_Traits::copy(this->_Bx()._Buf,
					_Right._Bx()._Buf, _Right._Mysize() + 1);
				this->_Getal().construct(::std:: addressof(_Right._Bx()._Ptr), _Ptr);
				}
		else
			if (_Right._Myres() < this->_BUF_SIZE)
				::std:: swap(this->_Bx()._Buf, _Right._Bx()._Buf);
			else
				{	
				pointer _Ptr = _Right._Bx()._Ptr;
				this->_Getal().destroy(::std:: addressof(_Right._Bx()._Ptr));
				_Traits::copy(_Right._Bx()._Buf,
					this->_Bx()._Buf, this->_Mysize() + 1);
				this->_Getal().construct(::std:: addressof(this->_Bx()._Ptr), _Ptr);
				}
		}

	void swap(_Myt& _Right)
		noexcept(_Alty::propagate_on_container_swap::value || _Alty::is_always_equal::value)
		{	
		if (this != &_Right)
			{	
			_Pocs(this->_Getal(), _Right._Getal());
			this->_Swap_all(_Right);
			_Swap_bx(_Right);
			::std:: swap(this->_Mysize(), _Right._Mysize());
			::std:: swap(this->_Myres(), _Right._Myres());
			}
		}

	size_type find(const _Myt& _Right, size_type _Off = 0) const noexcept
		{	
		return (find(_Right._Myptr(), _Off, _Right.size()));
		}

	size_type find(const _Elem *_Ptr,
		size_type _Off, size_type _Count) const
		{	
		;
		if (_Count == 0 && _Off <= this->_Mysize())
			return (_Off);	

		size_type _Nm;
		if (_Off < this->_Mysize() && _Count <= (_Nm = this->_Mysize() - _Off))
			{	
			const _Elem *_Uptr, *_Vptr;
			for (_Nm -= _Count - 1, _Vptr = this->_Myptr() + _Off;
				(_Uptr = _Traits::find(_Vptr, _Nm, *_Ptr)) != 0;
				_Nm -= _Uptr - _Vptr + 1, _Vptr = _Uptr + 1)
				if (_Traits::compare(_Uptr, _Ptr, _Count) == 0)
					return (_Uptr - this->_Myptr());	
			}

		return (npos);	
		}

	size_type find(const _Elem *_Ptr, size_type _Off = 0) const
		{	
		;
		return (find(_Ptr, _Off, _Traits::length(_Ptr)));
		}

	size_type find(_Elem _Ch, size_type _Off = 0) const
		{	
		return (find((const _Elem *)&_Ch, _Off, 1));
		}

	size_type rfind(const _Myt& _Right, size_type _Off = npos) const noexcept
		{	
		return (rfind(_Right._Myptr(), _Off, _Right.size()));
		}

	size_type rfind(const _Elem *_Ptr,
		size_type _Off, size_type _Count) const
		{	
		;
		if (_Count == 0)
			return (_Off < this->_Mysize() ? _Off
				: this->_Mysize());	
		if (_Count <= this->_Mysize())
			{	
			const _Elem *_Uptr = this->_Myptr() +
				(_Off < this->_Mysize() - _Count ? _Off
					: this->_Mysize() - _Count);
			for (; ; --_Uptr)
				if (_Traits::eq(*_Uptr, *_Ptr)
					&& _Traits::compare(_Uptr, _Ptr, _Count) == 0)
					return (_Uptr - this->_Myptr());	
				else if (_Uptr == this->_Myptr())
					break;	
			}

		return (npos);	
		}

	size_type rfind(const _Elem *_Ptr, size_type _Off = npos) const
		{	
		;
		return (rfind(_Ptr, _Off, _Traits::length(_Ptr)));
		}

	size_type rfind(_Elem _Ch, size_type _Off = npos) const
		{	
		return (rfind((const _Elem *)&_Ch, _Off, 1));
		}

	size_type find_first_of(const _Myt& _Right,
		size_type _Off = 0) const noexcept
		{	
		return (find_first_of(_Right._Myptr(), _Off, _Right.size()));
		}

	size_type find_first_of(const _Elem *_Ptr,
		size_type _Off, size_type _Count) const
		{	
		;
		if (0 < _Count && _Off < this->_Mysize())
			{	
			const _Elem *const _Vptr = this->_Myptr() + this->_Mysize();
			for (const _Elem *_Uptr = this->_Myptr() + _Off;
				_Uptr < _Vptr; ++_Uptr)
				if (_Traits::find(_Ptr, _Count, *_Uptr) != 0)
					return (_Uptr - this->_Myptr());	
			}

		return (npos);	
		}

	size_type find_first_of(const _Elem *_Ptr,
		size_type _Off = 0) const
		{	
		;
		return (find_first_of(_Ptr, _Off, _Traits::length(_Ptr)));
		}

	size_type find_first_of(_Elem _Ch,
		size_type _Off = 0) const
		{	
		return (find((const _Elem *)&_Ch, _Off, 1));
		}

	size_type find_last_of(const _Myt& _Right,
		size_type _Off = npos) const noexcept
		{	
		return (find_last_of(_Right._Myptr(), _Off, _Right.size()));
		}

	size_type find_last_of(const _Elem *_Ptr,
		size_type _Off, size_type _Count) const
		{	
		;
		if (0 < _Count && 0 < this->_Mysize())
			{	
			const _Elem *_Uptr = this->_Myptr()
				+ (_Off < this->_Mysize() ? _Off : this->_Mysize() - 1);
			for (; ; --_Uptr)
				if (_Traits::find(_Ptr, _Count, *_Uptr) != 0)
					return (_Uptr - this->_Myptr());	
				else if (_Uptr == this->_Myptr())
					break;	
			}

		return (npos);	
		}

	size_type find_last_of(const _Elem *_Ptr,
		size_type _Off = npos) const
		{	
		;
		return (find_last_of(_Ptr, _Off, _Traits::length(_Ptr)));
		}

	size_type find_last_of(_Elem _Ch,
		size_type _Off = npos) const
		{	
		return (rfind((const _Elem *)&_Ch, _Off, 1));
		}

	size_type find_first_not_of(const _Myt& _Right,
		size_type _Off = 0) const noexcept
		{	
		return (find_first_not_of(_Right._Myptr(), _Off,
			_Right.size()));
		}

	size_type find_first_not_of(const _Elem *_Ptr,
		size_type _Off, size_type _Count) const
		{	
		;
		if (_Off < this->_Mysize())
			{	
			const _Elem *const _Vptr = this->_Myptr() + this->_Mysize();
			for (const _Elem *_Uptr = this->_Myptr() + _Off;
				_Uptr < _Vptr; ++_Uptr)
				if (_Traits::find(_Ptr, _Count, *_Uptr) == 0)
					return (_Uptr - this->_Myptr());
			}
		return (npos);
		}

	size_type find_first_not_of(const _Elem *_Ptr,
		size_type _Off = 0) const
		{	
		;
		return (find_first_not_of(_Ptr, _Off, _Traits::length(_Ptr)));
		}

	size_type find_first_not_of(_Elem _Ch,
		size_type _Off = 0) const
		{	
		return (find_first_not_of((const _Elem *)&_Ch, _Off, 1));
		}

	size_type find_last_not_of(const _Myt& _Right,
		size_type _Off = npos) const noexcept
		{	
		return (find_last_not_of(_Right._Myptr(), _Off, _Right.size()));
		}

	size_type find_last_not_of(const _Elem *_Ptr,
		size_type _Off, size_type _Count) const
		{	
		;
		if (0 < this->_Mysize())
			{	
			const _Elem *_Uptr = this->_Myptr()
				+ (_Off < this->_Mysize() ? _Off : this->_Mysize() - 1);
			for (; ; --_Uptr)
				if (_Traits::find(_Ptr, _Count, *_Uptr) == 0)
					return (_Uptr - this->_Myptr());
				else if (_Uptr == this->_Myptr())
					break;
			}
		return (npos);
		}

	size_type find_last_not_of(const _Elem *_Ptr,
		size_type _Off = npos) const
		{	
		;
		return (find_last_not_of(_Ptr, _Off, _Traits::length(_Ptr)));
		}

	size_type find_last_not_of(_Elem _Ch,
		size_type _Off = npos) const
		{	
		return (find_last_not_of((const _Elem *)&_Ch, _Off, 1));
		}

	_Myt substr(size_type _Off = 0, size_type _Count = npos) const
		{	
		return (_Myt(*this, _Off, _Count, get_allocator()));
		}

	static int _Traits_compare(const _Elem * const _Left, const size_type _Left_size,
		const _Elem * const _Right, const size_type _Right_size)
		{	
		const size_type _Min_size = _Left_size < _Right_size ? _Left_size : _Right_size;
		const int _Ans = _Traits::compare(_Left, _Right, _Min_size);

		if (_Ans != 0)
			return (_Ans);

		if (_Left_size < _Right_size)
			return (-1);

		if (_Left_size > _Right_size)
			return (1);

		return (0);
		}

	size_type _Clamp_suffix_size(const size_type _Off, const size_type _Size) const
		{	
		const size_type _Max_effective_size = this->_Mysize() - _Off;
		if (_Size <= _Max_effective_size)
			return (_Size);
		else
			return (_Max_effective_size);
		}

	int compare(const _Myt& _Right) const noexcept
		{	
		return (_Traits_compare(this->_Myptr(), this->_Mysize(),
			_Right._Myptr(), _Right._Mysize()));
		}

	int compare(size_type _Off,
		size_type _N0, const _Myt& _Right) const
		{	
		_Check_offset(_Off);
		return (_Traits_compare(this->_Myptr() + _Off, _Clamp_suffix_size(_Off, _N0),
			_Right._Myptr(), _Right._Mysize()));
		}

	int compare(size_type _Off,
		size_type _N0, const _Myt& _Right,
		size_type _Roff, size_type _Count = npos) const
		{	
		_Check_offset(_Off);
		_Right._Check_offset(_Roff);
		return (_Traits_compare(this->_Myptr() + _Off, _Clamp_suffix_size(_Off, _N0),
			_Right._Myptr() + _Roff, _Right._Clamp_suffix_size(_Roff, _Count)));
		}

	int compare(const _Elem *_Ptr) const
		{	
		;
		return (_Traits_compare(this->_Myptr(), this->_Mysize(),
			_Ptr, _Traits::length(_Ptr)));
		}

	int compare(size_type _Off, size_type _N0, const _Elem *_Ptr) const
		{	
		;
		_Check_offset(_Off);
		return (_Traits_compare(this->_Myptr() + _Off, _Clamp_suffix_size(_Off, _N0),
			_Ptr, _Traits::length(_Ptr)));
		}

	int compare(size_type _Off,
		size_type _N0, const _Elem *_Ptr, size_type _Count) const
		{	
		;
		_Check_offset(_Off);
		return (_Traits_compare(this->_Myptr() + _Off, _Clamp_suffix_size(_Off, _N0),
			_Ptr, _Count));
		}

	allocator_type get_allocator() const noexcept
		{	
		allocator_type _Ret(this->_Getal());
		return (_Ret);
		}

	void _Chassign(size_type _Off, size_type _Count, _Elem _Ch)
		{	
		if (_Count == 1)
			_Traits::assign(*(this->_Myptr() + _Off), _Ch);
		else
			_Traits::assign(this->_Myptr() + _Off, _Count, _Ch);
		}

	void _Copy(size_type _Newsize, size_type _Oldlen)
		{	
		size_type _Newres = _Newsize | this->_ALLOC_MASK;
		if (max_size() < _Newres)
			_Newres = _Newsize;	
		else if (this->_Myres() / 2 <= _Newres / 3)
			;
		else if (this->_Myres() <= max_size() - this->_Myres() / 2)
			_Newres = this->_Myres()
				+ this->_Myres() / 2;	
		else
			_Newres = max_size();	

		pointer _Ptr;
		try {
			_Ptr = this->_Getal().allocate(_Newres + 1);
		} catch (...) {
			_Newres = _Newsize;	
			try {
				_Ptr = this->_Getal().allocate(_Newres + 1);
			} catch (...) {
			_Tidy(true);	
			throw;
			}
		}

		if (0 < _Oldlen)
			_Traits::copy(_Unfancy(_Ptr), this->_Myptr(),
				_Oldlen);	
		_Tidy(true);
		this->_Getal().construct(::std:: addressof(this->_Bx()._Ptr), _Ptr);
		this->_Myres() = _Newres;
		_Eos(_Oldlen);
		}

	void _Eos(size_type _Newsize)
		{	
		auto& _Dx = this->_Get_data();
		_Traits::assign(_Dx._Myptr()[_Dx._Mysize = _Newsize], _Elem());
		}

	bool _Grow(size_type _Newsize,
		bool _Trim = false)
		{	
		if (max_size() < _Newsize)
			_Xlen();	
		if (this->_Myres() < _Newsize)
			_Copy(_Newsize, this->_Mysize());	
		else if (_Trim && _Newsize < this->_BUF_SIZE)
			_Tidy(true,	
				_Newsize < this->_Mysize() ? _Newsize : this->_Mysize());
		else if (_Newsize == 0)
			_Eos(0);	
		return (0 < _Newsize);	
		}

	bool _Inside(const _Elem *_Ptr)
		{	
		if (_Ptr == nullptr_t{} || _Ptr < this->_Myptr()
			|| this->_Myptr() + this->_Mysize() <= _Ptr)
			return (false);	
		else
			return (true);
		}

	void _Tidy(bool _Built = false,
		size_type _Newsize = 0)
		{	
		if (!_Built)
			;
		else if (this->_BUF_SIZE <= this->_Myres())
			{	
			pointer _Ptr = this->_Bx()._Ptr;
			this->_Getal().destroy(::std:: addressof(this->_Bx()._Ptr));
			if (0 < _Newsize)
				_Traits::copy(this->_Bx()._Buf,
					_Unfancy(_Ptr), _Newsize);
			this->_Getal().deallocate(_Ptr, this->_Myres() + 1);
			}
		this->_Myres() = this->_BUF_SIZE - 1;
		_Eos(_Newsize);
		}

	[[noreturn]] void _Xlen() const
		{	
		_Xlength_error("string too long");
		}

	void _Check_offset(const size_type _Off) const
		{	
		if (this->_Mysize() < _Off)
			_Xran();
		}

	void _Check_offset_exclusive(const size_type _Off) const
		{	
		if (this->_Mysize() <= _Off)
			_Xran();
		}

	[[noreturn]] void _Xran() const
		{	
		_Xout_of_range("invalid string position");
		}
	};

		
template<class _Elem,
	class _Traits,
	class _Alloc>
	 const typename basic_string<_Elem, _Traits, _Alloc>::size_type
		basic_string<_Elem, _Traits, _Alloc>::npos =
			(typename basic_string<_Elem, _Traits, _Alloc>::size_type)(-1);

		

template<class _Elem,
	class _Traits,
	class _Alloc> inline
	void swap(basic_string<_Elem, _Traits, _Alloc>& _Left,
		basic_string<_Elem, _Traits, _Alloc>& _Right)
			noexcept(noexcept(_Left.swap(_Right)))
	{	
	_Left.swap(_Right);
	}

template<class _Elem,
	class _Traits,
	class _Alloc> inline
	basic_string<_Elem, _Traits, _Alloc> operator+(
		const basic_string<_Elem, _Traits, _Alloc>& _Left,
		const basic_string<_Elem, _Traits, _Alloc>& _Right)
	{	
	basic_string<_Elem, _Traits, _Alloc> _Ans;
	_Ans.reserve(_Left.size() + _Right.size());
	_Ans += _Left;
	_Ans += _Right;
	return (_Ans);
	}

template<class _Elem,
	class _Traits,
	class _Alloc> inline
	basic_string<_Elem, _Traits, _Alloc> operator+(
		const _Elem *_Left,
		const basic_string<_Elem, _Traits, _Alloc>& _Right)
	{	
	basic_string<_Elem, _Traits, _Alloc> _Ans;
	_Ans.reserve(_Traits::length(_Left) + _Right.size());
	_Ans += _Left;
	_Ans += _Right;
	return (_Ans);
	}

template<class _Elem,
	class _Traits,
	class _Alloc> inline
	basic_string<_Elem, _Traits, _Alloc> operator+(
		const _Elem _Left,
		const basic_string<_Elem, _Traits, _Alloc>& _Right)
	{	
	basic_string<_Elem, _Traits, _Alloc> _Ans;
	_Ans.reserve(1 + _Right.size());
	_Ans += _Left;
	_Ans += _Right;
	return (_Ans);
	}

template<class _Elem,
	class _Traits,
	class _Alloc> inline
	basic_string<_Elem, _Traits, _Alloc> operator+(
		const basic_string<_Elem, _Traits, _Alloc>& _Left,
		const _Elem *_Right)
	{	
	basic_string<_Elem, _Traits, _Alloc> _Ans;
	_Ans.reserve(_Left.size() + _Traits::length(_Right));
	_Ans += _Left;
	_Ans += _Right;
	return (_Ans);
	}

template<class _Elem,
	class _Traits,
	class _Alloc> inline
	basic_string<_Elem, _Traits, _Alloc> operator+(
		const basic_string<_Elem, _Traits, _Alloc>& _Left,
		const _Elem _Right)
	{	
	basic_string<_Elem, _Traits, _Alloc> _Ans;
	_Ans.reserve(_Left.size() + 1);
	_Ans += _Left;
	_Ans += _Right;
	return (_Ans);
	}

template<class _Elem,
	class _Traits,
	class _Alloc> inline
	basic_string<_Elem, _Traits, _Alloc> operator+(
		const basic_string<_Elem, _Traits, _Alloc>& _Left,
		basic_string<_Elem, _Traits, _Alloc>&& _Right)
	{	
	return (::std:: move(_Right.insert(0, _Left)));
	}

template<class _Elem,
	class _Traits,
	class _Alloc> inline
	basic_string<_Elem, _Traits, _Alloc> operator+(
		basic_string<_Elem, _Traits, _Alloc>&& _Left,
		const basic_string<_Elem, _Traits, _Alloc>& _Right)
	{	
	return (::std:: move(_Left.append(_Right)));
	}

template<class _Elem,
	class _Traits,
	class _Alloc> inline
	basic_string<_Elem, _Traits, _Alloc> operator+(
		basic_string<_Elem, _Traits, _Alloc>&& _Left,
		basic_string<_Elem, _Traits, _Alloc>&& _Right)
	{	
	if (_Right.size() <= _Left.capacity() - _Left.size()
		|| _Right.capacity() - _Right.size() < _Left.size())
		return (::std:: move(_Left.append(_Right)));
	else
		return (::std:: move(_Right.insert(0, _Left)));
	}

template<class _Elem,
	class _Traits,
	class _Alloc> inline
	basic_string<_Elem, _Traits, _Alloc> operator+(
		const _Elem *_Left,
		basic_string<_Elem, _Traits, _Alloc>&& _Right)
	{	
	return (::std:: move(_Right.insert(0, _Left)));
	}

template<class _Elem,
	class _Traits,
	class _Alloc> inline
	basic_string<_Elem, _Traits, _Alloc> operator+(
		const _Elem _Left,
		basic_string<_Elem, _Traits, _Alloc>&& _Right)
	{	
	typedef typename basic_string<_Elem, _Traits, _Alloc>::size_type
		size_type;
	return (::std:: move(_Right.insert((size_type)0, (size_type)1, _Left)));
	}

template<class _Elem,
	class _Traits,
	class _Alloc> inline
	basic_string<_Elem, _Traits, _Alloc> operator+(
		basic_string<_Elem, _Traits, _Alloc>&& _Left,
		const _Elem *_Right)
	{	
	return (::std:: move(_Left.append(_Right)));
	}

template<class _Elem,
	class _Traits,
	class _Alloc> inline
	basic_string<_Elem, _Traits, _Alloc> operator+(
		basic_string<_Elem, _Traits, _Alloc>&& _Left,
		const _Elem _Right)
	{	
	return (::std:: move(_Left.append(1, _Right)));
	}

template<class _Elem,
	class _Traits,
	class _Alloc> inline
	bool operator==(
		const basic_string<_Elem, _Traits, _Alloc>& _Left,
		const basic_string<_Elem, _Traits, _Alloc>& _Right) noexcept
	{	
	return (_Left.compare(_Right) == 0);
	}

template<class _Elem,
	class _Traits,
	class _Alloc> inline
	bool operator==(
		const _Elem * _Left,
		const basic_string<_Elem, _Traits, _Alloc>& _Right)
	{	
	return (_Right.compare(_Left) == 0);
	}

template<class _Elem,
	class _Traits,
	class _Alloc> inline
	bool operator==(
		const basic_string<_Elem, _Traits, _Alloc>& _Left,
		const _Elem *_Right)
	{	
	return (_Left.compare(_Right) == 0);
	}

template<class _Elem,
	class _Traits,
	class _Alloc> inline
	bool operator!=(
		const basic_string<_Elem, _Traits, _Alloc>& _Left,
		const basic_string<_Elem, _Traits, _Alloc>& _Right) noexcept
	{	
	return (!(_Left == _Right));
	}

template<class _Elem,
	class _Traits,
	class _Alloc> inline
	bool operator!=(
		const _Elem *_Left,
		const basic_string<_Elem, _Traits, _Alloc>& _Right)
	{	
	return (!(_Left == _Right));
	}

template<class _Elem,
	class _Traits,
	class _Alloc> inline
	bool operator!=(
		const basic_string<_Elem, _Traits, _Alloc>& _Left,
		const _Elem *_Right)
	{	
	return (!(_Left == _Right));
	}

template<class _Elem,
	class _Traits,
	class _Alloc> inline
	bool operator<(
		const basic_string<_Elem, _Traits, _Alloc>& _Left,
		const basic_string<_Elem, _Traits, _Alloc>& _Right) noexcept
	{	
	return (_Left.compare(_Right) < 0);
	}

template<class _Elem,
	class _Traits,
	class _Alloc> inline
	bool operator<(
		const _Elem * _Left,
		const basic_string<_Elem, _Traits, _Alloc>& _Right)
	{	
	return (_Right.compare(_Left) > 0);
	}

template<class _Elem,
	class _Traits,
	class _Alloc> inline
	bool operator<(
		const basic_string<_Elem, _Traits, _Alloc>& _Left,
		const _Elem *_Right)
	{	
	return (_Left.compare(_Right) < 0);
	}

template<class _Elem,
	class _Traits,
	class _Alloc> inline
	bool operator>(
		const basic_string<_Elem, _Traits, _Alloc>& _Left,
		const basic_string<_Elem, _Traits, _Alloc>& _Right) noexcept
	{	
	return (_Right < _Left);
	}

template<class _Elem,
	class _Traits,
	class _Alloc> inline
	bool operator>(
		const _Elem * _Left,
		const basic_string<_Elem, _Traits, _Alloc>& _Right)
	{	
	return (_Right < _Left);
	}

template<class _Elem,
	class _Traits,
	class _Alloc> inline
	bool operator>(
		const basic_string<_Elem, _Traits, _Alloc>& _Left,
		const _Elem *_Right)
	{	
	return (_Right < _Left);
	}

template<class _Elem,
	class _Traits,
	class _Alloc> inline
	bool operator<=(
		const basic_string<_Elem, _Traits, _Alloc>& _Left,
		const basic_string<_Elem, _Traits, _Alloc>& _Right) noexcept
	{	
	return (!(_Right < _Left));
	}

template<class _Elem,
	class _Traits,
	class _Alloc> inline
	bool operator<=(
		const _Elem * _Left,
		const basic_string<_Elem, _Traits, _Alloc>& _Right)
	{	
	return (!(_Right < _Left));
	}

template<class _Elem,
	class _Traits,
	class _Alloc> inline
	bool operator<=(
		const basic_string<_Elem, _Traits, _Alloc>& _Left,
		const _Elem *_Right)
	{	
	return (!(_Right < _Left));
	}

template<class _Elem,
	class _Traits,
	class _Alloc> inline
	bool operator>=(
		const basic_string<_Elem, _Traits, _Alloc>& _Left,
		const basic_string<_Elem, _Traits, _Alloc>& _Right) noexcept
	{	
	return (!(_Left < _Right));
	}

template<class _Elem,
	class _Traits,
	class _Alloc> inline
	bool operator>=(
		const _Elem * _Left,
		const basic_string<_Elem, _Traits, _Alloc>& _Right)
	{	
	return (!(_Left < _Right));
	}

template<class _Elem,
	class _Traits,
	class _Alloc> inline
	bool operator>=(
		const basic_string<_Elem, _Traits, _Alloc>& _Left,
		const _Elem *_Right)
	{	
	return (!(_Left < _Right));
	}

typedef basic_string<char, char_traits<char>, allocator<char> >
	string;
typedef basic_string<wchar_t, char_traits<wchar_t>, allocator<wchar_t> >
	wstring;

	
template<class _Elem,
	class _Traits,
	class _Alloc>
	struct hash<basic_string<_Elem, _Traits, _Alloc> >
	{	
	typedef basic_string<_Elem, _Traits, _Alloc> argument_type;
	typedef size_t result_type;

	size_t operator()(const argument_type& _Keyval) const
		{	
		return (_Hash_seq((const unsigned char *)_Keyval.c_str(),
			_Keyval.size() * sizeof (_Elem)));
		}
	};

typedef basic_string<char16_t, char_traits<char16_t>, allocator<char16_t> >
	u16string;
typedef basic_string<char32_t, char_traits<char32_t>, allocator<char32_t> >
	u32string;
}

 
 #pragma warning(pop)
 #pragma pack(pop)










 #pragma pack(push,8)
 #pragma warning(push,3)
 
 
namespace std {
		
class logic_error
	: public ::std:: exception
	{	
public:
	typedef ::std:: exception _Mybase;

	explicit logic_error(const string& _Message)
		: _Mybase(_Message.c_str())
		{	
		}

	explicit logic_error(const char *_Message)
		: _Mybase(_Message)
		{	
		}

 

 






	};

		
class domain_error
	: public logic_error
	{	
public:
	typedef logic_error _Mybase;

	explicit domain_error(const string& _Message)
		: _Mybase(_Message.c_str())
		{	
		}

	explicit domain_error(const char *_Message)
		: _Mybase(_Message)
		{	
		}

 

 






	};

		
class invalid_argument
	: public logic_error
	{	
public:
	typedef logic_error _Mybase;

	explicit invalid_argument(const string& _Message)
		: _Mybase(_Message.c_str())
		{	
		}

	explicit invalid_argument(const char *_Message)
		: _Mybase(_Message)
		{	
		}

 

 






	};

		
class length_error
	: public logic_error
	{	
public:
	typedef logic_error _Mybase;

	explicit length_error(const string& _Message)
		: _Mybase(_Message.c_str())
		{	
		}

	explicit length_error(const char *_Message)
		: _Mybase(_Message)
		{	
		}

 

 






	};

		
class out_of_range
	: public logic_error
	{	
public:
	typedef logic_error _Mybase;

	explicit out_of_range(const string& _Message)
		: _Mybase(_Message.c_str())
		{	
		}

	explicit out_of_range(const char *_Message)
		: _Mybase(_Message)
		{	
		}

 

 






	};

		
class runtime_error
	: public ::std:: exception
	{	
public:
	typedef ::std:: exception _Mybase;

	explicit runtime_error(const string& _Message)
		: _Mybase(_Message.c_str())
		{	
		}

	explicit runtime_error(const char *_Message)
		: _Mybase(_Message)
		{	
		}

 

 






	};

		
class overflow_error
	: public runtime_error
	{	
public:
	typedef runtime_error _Mybase;

	explicit overflow_error(const string& _Message)
		: _Mybase(_Message.c_str())
		{	
		}

	explicit overflow_error(const char *_Message)
		: _Mybase(_Message)
		{	
		}

 

 






	};

		
class underflow_error
	: public runtime_error
	{	
public:
	typedef runtime_error _Mybase;

	explicit underflow_error(const string& _Message)
		: _Mybase(_Message.c_str())
		{	
		}

	explicit underflow_error(const char *_Message)
		: _Mybase(_Message)
		{	
		}

 

 






	};

		
class range_error
	: public runtime_error
	{	
public:
	typedef runtime_error _Mybase;

	explicit range_error(const string& _Message)
		: _Mybase(_Message.c_str())
		{	
		}

	explicit range_error(const char *_Message)
		: _Mybase(_Message)
		{	
		}

 

 






	};
}
 
 #pragma warning(pop)
 #pragma pack(pop)









 #pragma pack(push,8)
 #pragma warning(push,3)
 
 

 #pragma warning(disable: 4127)
 #pragma warning(disable: 4244)

namespace std {
 

		
template<class _Myvec>
	class _Vector_const_iterator
		: public _Iterator012<random_access_iterator_tag,
			typename _Myvec::value_type,
			typename _Myvec::difference_type,
			typename _Myvec::const_pointer,
			typename _Myvec::const_reference,
			_Iterator_base>
	{	
public:
	typedef _Vector_const_iterator<_Myvec> _Myiter;
	typedef random_access_iterator_tag iterator_category;

	typedef typename _Myvec::value_type value_type;
	typedef typename _Myvec::difference_type difference_type;
	typedef typename _Myvec::const_pointer pointer;
	typedef typename _Myvec::const_reference reference;
	typedef typename _Myvec::pointer _Tptr;

	_Vector_const_iterator()
		: _Ptr()
		{	
		}

	_Vector_const_iterator(_Tptr _Parg, const _Container_base *_Pvector)
		: _Ptr(_Parg)
		{	
		this->_Adopt(_Pvector);
		}

	typedef pointer _Unchecked_type;

	_Myiter& _Rechecked(_Unchecked_type _Right)
		{	
		_Ptr = _Const_cast(_Right);
		return (*this);
		}

	_Unchecked_type _Unchecked() const
		{	
		return (_Ptr);
		}

	reference operator*() const
		{	
 

















		;

		return (*_Ptr);
		}

	pointer operator->() const
		{	
		return (pointer_traits<pointer>::pointer_to(**this));
		}

	_Myiter& operator++()
		{	
 
















		++_Ptr;
		return (*this);
		}

	_Myiter operator++(int)
		{	
		_Myiter _Tmp = *this;
		++*this;
		return (_Tmp);
		}

	_Myiter& operator--()
		{	
 
















		--_Ptr;
		return (*this);
		}

	_Myiter operator--(int)
		{	
		_Myiter _Tmp = *this;
		--*this;
		return (_Tmp);
		}

	_Myiter& operator+=(difference_type _Off)
		{	
 






















		_Ptr += _Off;
		return (*this);
		}

	_Myiter operator+(difference_type _Off) const
		{	
		_Myiter _Tmp = *this;
		return (_Tmp += _Off);
		}

	_Myiter& operator-=(difference_type _Off)
		{	
		return (*this += -_Off);
		}

	_Myiter operator-(difference_type _Off) const
		{	
		_Myiter _Tmp = *this;
		return (_Tmp -= _Off);
		}

	difference_type operator-(const _Myiter& _Right) const
		{	
		_Compat(_Right);
		return (_Ptr - _Right._Ptr);
		}

	reference operator[](difference_type _Off) const
		{	
		return (*(*this + _Off));
		}

	bool operator==(const _Myiter& _Right) const
		{	
		_Compat(_Right);
		return (_Ptr == _Right._Ptr);
		}

	bool operator!=(const _Myiter& _Right) const
		{	
		return (!(*this == _Right));
		}

	bool operator<(const _Myiter& _Right) const
		{	
		_Compat(_Right);
		return (_Ptr < _Right._Ptr);
		}

	bool operator>(const _Myiter& _Right) const
		{	
		return (_Right < *this);
		}

	bool operator<=(const _Myiter& _Right) const
		{	
		return (!(_Right < *this));
		}

	bool operator>=(const _Myiter& _Right) const
		{	
		return (!(*this < _Right));
		}

 
















	void _Compat(const _Myiter&) const
		{	
		}
 

	_Tptr _Ptr;	
	};

template<class _Myvec> inline
	typename _Vector_const_iterator<_Myvec>::_Unchecked_type
		_Unchecked(_Vector_const_iterator<_Myvec> _Iter)
	{	
	return (_Iter._Unchecked());
	}

template<class _Myvec> inline
	_Vector_const_iterator<_Myvec>&
		_Rechecked(_Vector_const_iterator<_Myvec>& _Iter,
			typename _Vector_const_iterator<_Myvec>
				::_Unchecked_type _Right)
	{	
	return (_Iter._Rechecked(_Right));
	}

template<class _Myvec> inline
	_Vector_const_iterator<_Myvec> operator+(
		typename _Vector_const_iterator<_Myvec>::difference_type _Off,
		_Vector_const_iterator<_Myvec> _Next)
	{	
	return (_Next += _Off);
	}

		
template<class _Myvec>
	class _Vector_iterator
		: public _Vector_const_iterator<_Myvec>
	{	
public:
	typedef _Vector_iterator<_Myvec> _Myiter;
	typedef _Vector_const_iterator<_Myvec> _Mybase;
	typedef random_access_iterator_tag iterator_category;

	typedef typename _Myvec::value_type value_type;
	typedef typename _Myvec::difference_type difference_type;
	typedef typename _Myvec::pointer pointer;
	typedef typename _Myvec::reference reference;

	_Vector_iterator()
		{	
		}

	_Vector_iterator(pointer _Parg, const _Container_base *_Pvector)
		: _Mybase(_Parg, _Pvector)
		{	
		}

	typedef pointer _Unchecked_type;

	_Myiter& _Rechecked(_Unchecked_type _Right)
		{	
		this->_Ptr = _Right;
		return (*this);
		}

	_Unchecked_type _Unchecked() const
		{	
		return (this->_Ptr);
		}

	reference operator*() const
		{	
		return ((reference)**(_Mybase *)this);
		}

	pointer operator->() const
		{	
		return (pointer_traits<pointer>::pointer_to(**this));
		}

	_Myiter& operator++()
		{	
		++*(_Mybase *)this;
		return (*this);
		}

	_Myiter operator++(int)
		{	
		_Myiter _Tmp = *this;
		++*this;
		return (_Tmp);
		}

	_Myiter& operator--()
		{	
		--*(_Mybase *)this;
		return (*this);
		}

	_Myiter operator--(int)
		{	
		_Myiter _Tmp = *this;
		--*this;
		return (_Tmp);
		}

	_Myiter& operator+=(difference_type _Off)
		{	
		*(_Mybase *)this += _Off;
		return (*this);
		}

	_Myiter operator+(difference_type _Off) const
		{	
		_Myiter _Tmp = *this;
		return (_Tmp += _Off);
		}

	_Myiter& operator-=(difference_type _Off)
		{	
		return (*this += -_Off);
		}

	_Myiter operator-(difference_type _Off) const
		{	
		_Myiter _Tmp = *this;
		return (_Tmp -= _Off);
		}

	difference_type operator-(const _Mybase& _Right) const
		{	
		return (*(_Mybase *)this - _Right);
		}

	reference operator[](difference_type _Off) const
		{	
		return (*(*this + _Off));
		}
	};

template<class _Myvec> inline
	typename _Vector_iterator<_Myvec>::_Unchecked_type
		_Unchecked(_Vector_iterator<_Myvec> _Iter)
	{	
	return (_Iter._Unchecked());
	}

template<class _Myvec> inline
	_Vector_iterator<_Myvec>&
		_Rechecked(_Vector_iterator<_Myvec>& _Iter,
			typename _Vector_iterator<_Myvec>
				::_Unchecked_type _Right)
	{	
	return (_Iter._Rechecked(_Right));
	}

template<class _Myvec> inline
	_Vector_iterator<_Myvec> operator+(
		typename _Vector_iterator<_Myvec>::difference_type _Off,
		_Vector_iterator<_Myvec> _Next)
	{	
	return (_Next += _Off);
	}

		
template<class _Value_type,
	class _Size_type,
	class _Difference_type,
	class _Pointer,
	class _Const_pointer,
	class _Reference,
	class _Const_reference>
	struct _Vec_iter_types
	{	
	typedef _Value_type value_type;
	typedef _Size_type size_type;
	typedef _Difference_type difference_type;
	typedef _Pointer pointer;
	typedef _Const_pointer const_pointer;
	typedef _Reference reference;
	typedef _Const_reference const_reference;
	};

template<class _Ty,
	class _Alloc0>
	struct _Vec_base_types
	{	
	typedef _Alloc0 _Alloc;
	typedef _Vec_base_types<_Ty, _Alloc> _Myt;

	typedef _Wrap_alloc<_Alloc> _Alty0;
	typedef typename _Alty0::template rebind<_Ty>::other _Alty;


	typedef typename _If<_Is_simple_alloc<_Alty>::value,
		_Simple_types<typename _Alty::value_type>,
		_Vec_iter_types<typename _Alty::value_type,
			typename _Alty::size_type,
			typename _Alty::difference_type,
			typename _Alty::pointer,
			typename _Alty::const_pointer,
			typename _Alty::reference,
			typename _Alty::const_reference> >::type
		_Val_types;
	};

		
template<class _Val_types>
	class _Vector_val
		: public _Container_base
	{	
public:
	typedef _Vector_val<_Val_types> _Myt;

	typedef typename _Val_types::value_type value_type;
	typedef typename _Val_types::size_type size_type;
	typedef typename _Val_types::difference_type difference_type;
	typedef typename _Val_types::pointer pointer;
	typedef typename _Val_types::const_pointer const_pointer;
	typedef typename _Val_types::reference reference;
	typedef typename _Val_types::const_reference const_reference;

	typedef _Vector_iterator<_Myt> iterator;
	typedef _Vector_const_iterator<_Myt> const_iterator;

	_Vector_val()
		: _Myfirst(),
		_Mylast(),
		_Myend()
		{	
		}

	pointer _Myfirst;	
	pointer _Mylast;	
	pointer _Myend;	
	};

		
template<class _Alloc_types>
	class _Vector_alloc
	{	
public:
	typedef _Vector_alloc<_Alloc_types> _Myt;
	typedef typename _Alloc_types::_Alloc _Alloc;
	typedef typename _Alloc_types::_Alty _Alty;
	typedef typename _Alloc_types::_Val_types _Val_types;

	typedef typename _Val_types::value_type value_type;
	typedef typename _Val_types::size_type size_type;
	typedef typename _Val_types::difference_type difference_type;
	typedef typename _Val_types::pointer pointer;
	typedef typename _Val_types::const_pointer const_pointer;
	typedef typename _Val_types::reference reference;
	typedef typename _Val_types::const_reference const_reference;

	typedef _Vector_iterator<_Vector_val<_Val_types> > iterator;
	typedef _Vector_const_iterator<_Vector_val<_Val_types> > const_iterator;

 
	_Vector_alloc()
		: _Mypair(_Zero_then_variadic_args_t())
		{	
		}

	template<class _Any_alloc,
		class = enable_if_t<!is_same<decay_t<_Any_alloc>, _Myt>::value> >
		_Vector_alloc(_Any_alloc&& _Al)
		: _Mypair(_One_then_variadic_args_t(),
			::std:: forward<_Any_alloc>(_Al))
		{	
		}

	void _Copy_alloc(const _Alty& _Al)
		{	
		_Pocca(_Getal(), _Al);
		}

	void _Move_alloc(_Alty& _Al)
		{	
		_Pocma(_Getal(), _Al);
		}

 





































































	void _Orphan_all()
		{	
		_Get_data()._Orphan_all();
		}

	void _Swap_all(_Myt& _Right)
		{	
		_Get_data()._Swap_all(_Right._Get_data());
		}

	_Alty& _Getal() noexcept
		{	
		return (_Mypair._Get_first());
		}

	const _Alty& _Getal() const noexcept
		{	
		return (_Mypair._Get_first());
		}

	_Vector_val<_Val_types>& _Get_data() noexcept
		{	
		return (_Mypair._Get_second());
		}

	const _Vector_val<_Val_types>& _Get_data() const noexcept
		{	
		return (_Mypair._Get_second());
		}

	pointer& _Myfirst() noexcept
		{	
		return (_Get_data()._Myfirst);
		}

	const pointer& _Myfirst() const noexcept
		{	
		return (_Get_data()._Myfirst);
		}

	pointer& _Mylast() noexcept
		{	
		return (_Get_data()._Mylast);
		}

	const pointer& _Mylast() const noexcept
		{	
		return (_Get_data()._Mylast);
		}

	pointer& _Myend() noexcept
		{	
		return (_Get_data()._Myend);
		}

	const pointer& _Myend() const noexcept
		{	
		return (_Get_data()._Myend);
		}

private:
	_Compressed_pair<_Alty, _Vector_val<_Val_types> > _Mypair;
	};

		
template<class _Ty,
	class _Alloc = allocator<_Ty> >
	class vector
		: public _Vector_alloc<_Vec_base_types<_Ty, _Alloc> >
	{	
public:
	typedef vector<_Ty, _Alloc> _Myt;
	typedef _Vector_alloc<_Vec_base_types<_Ty, _Alloc> > _Mybase;
	typedef _Alloc allocator_type;

	typedef typename _Mybase::_Alty _Alty;

	typedef typename _Mybase::value_type value_type;
	typedef typename _Mybase::size_type size_type;
	typedef typename _Mybase::difference_type difference_type;
	typedef typename _Mybase::pointer pointer;
	typedef typename _Mybase::const_pointer const_pointer;
	typedef typename _Mybase::reference reference;
	typedef typename _Mybase::const_reference const_reference;

 
 

	typedef typename _Mybase::iterator iterator;
	typedef typename _Mybase::const_iterator const_iterator;

	typedef ::std:: reverse_iterator<iterator> reverse_iterator;
	typedef ::std:: reverse_iterator<const_iterator> const_reverse_iterator;

	vector() noexcept(is_nothrow_default_constructible<_Alloc>::value)
		: _Mybase()
		{	
		}

	explicit vector(const _Alloc& _Al) noexcept
		: _Mybase(_Al)
		{	
		}

	explicit vector(size_type _Count)
		: _Mybase()
		{	
		if (_Buy(_Count))
			{	
			try {
			_Uninitialized_default_fill_n(this->_Myfirst(), _Count,
				this->_Getal());
			this->_Mylast() += _Count;
			} catch (...) {
			_Tidy();
			throw;
			}
			}
		}

	vector(size_type _Count, const value_type& _Val)
		: _Mybase()
		{	
		_Construct_n(_Count, ::std:: addressof(_Val));
		}

	vector(size_type _Count, const value_type& _Val, const _Alloc& _Al)
		: _Mybase(_Al)
		{	
		_Construct_n(_Count, ::std:: addressof(_Val));
		}

	vector(const _Myt& _Right)

		: _Mybase(_Right._Getal().select_on_container_copy_construction())


		{	
		if (_Buy(_Right.size()))
			try {
			this->_Mylast() = _Ucopy(_Right.begin(), _Right.end(),
				this->_Myfirst());
			} catch (...) {
			_Tidy();
			throw;
			}
		}

	vector(const _Myt& _Right, const _Alloc& _Al)
		: _Mybase(_Al)
		{	
		if (_Buy(_Right.size()))
			try {
			this->_Mylast() = _Ucopy(_Right.begin(), _Right.end(),
				this->_Myfirst());
			} catch (...) {
			_Tidy();
			throw;
			}
		}

	template<class _Iter,
		class = typename enable_if<_Is_iterator<_Iter>::value,
			void>::type>
		vector(_Iter _First, _Iter _Last)
		: _Mybase()
		{	
		_Construct(_First, _Last);
		}

	template<class _Iter,
		class = typename enable_if<_Is_iterator<_Iter>::value,
			void>::type>
		vector(_Iter _First, _Iter _Last, const _Alloc& _Al)
		: _Mybase(_Al)
		{	
		_Construct(_First, _Last);
		}

	template<class _Iter>
		void _Construct(_Iter _First, _Iter _Last)
		{	
		_Construct(_First, _Last, _Iter_cat_t<_Iter>());
		}

	template<class _Iter>
		void _Construct(_Iter _First, _Iter _Last,
			input_iterator_tag)
		{	
		try {

		for (; _First != _Last; ++_First)
			emplace_back(*_First);

		} catch (...) {
		_Tidy();
		throw;
		}
		}

	template<class _Iter>
		void _Construct(_Iter _First, _Iter _Last,
			forward_iterator_tag)
		{	
		if (_Buy(::std:: distance(_First, _Last)))
			{	
			try {
			this->_Mylast() = _Ucopy(_First, _Last, this->_Myfirst());
			} catch (...) {
			_Tidy();
			throw;
			}
			}
		}

	void _Construct_n(size_type _Count, const value_type *_Pval)
		{	
		if (_Buy(_Count))
			{	
			try {
			this->_Mylast() = _Ufill(this->_Myfirst(), _Count, _Pval);
			} catch (...) {
			_Tidy();
			throw;
			}
			}
		}

	vector(_Myt&& _Right) noexcept
		: _Mybase(::std:: move(_Right._Getal()))
		{	
		_Assign_rv(::std:: forward<_Myt>(_Right), true_type());
		}

	vector(_Myt&& _Right, const _Alloc& _Al)
		: _Mybase(_Al)
		{	
		_Assign_rv(::std:: forward<_Myt>(_Right));
		}

	_Myt& operator=(_Myt&& _Right)
		noexcept(_Alty::propagate_on_container_move_assignment::value || _Alty::is_always_equal::value)
		{	
		if (this != &_Right)
			{	
			_Tidy();
			if (_Alty::propagate_on_container_move_assignment::value
				&& this->_Getal() != _Right._Getal())
				this->_Move_alloc(_Right._Getal());

			_Assign_rv(::std:: forward<_Myt>(_Right));
			}
		return (*this);
		}

	void _Assign_rv(_Myt&& _Right, true_type)
		{	
		this->_Swap_all((_Myt&)_Right);
		this->_Myfirst() = _Right._Myfirst();
		this->_Mylast() = _Right._Mylast();
		this->_Myend() = _Right._Myend();

		_Right._Myfirst() = pointer();
		_Right._Mylast() = pointer();
		_Right._Myend() = pointer();
		}

	void _Assign_rv(_Myt&& _Right, false_type)
		{	
		if (get_allocator() == _Right.get_allocator())
			_Assign_rv(::std:: forward<_Myt>(_Right), true_type());
		else
			_Construct(::std:: make_move_iterator(_Right.begin()),
				::std:: make_move_iterator(_Right.end()));
		}

	void _Assign_rv(_Myt&& _Right)
		{	
		_Assign_rv(::std:: forward<_Myt>(_Right),
			typename _Alty::propagate_on_container_move_assignment());
		}


	void push_back(value_type&& _Val)
		{	
		if (_Inside(::std:: addressof(_Val)))
			{	
			size_type _Idx = ::std:: addressof(_Val) - _Unfancy(this->_Myfirst());
			if (this->_Mylast() == this->_Myend())
				_Reserve(1);
			_Orphan_range(this->_Mylast(), this->_Mylast());
			this->_Getal().construct(_Unfancy(this->_Mylast()),
				::std:: forward<value_type>(this->_Myfirst()[_Idx]));
			++this->_Mylast();
			}
		else
			{	
			if (this->_Mylast() == this->_Myend())
				_Reserve(1);
			_Orphan_range(this->_Mylast(), this->_Mylast());
			this->_Getal().construct(_Unfancy(this->_Mylast()),
				::std:: forward<value_type>(_Val));
			++this->_Mylast();
			}
		}

	iterator insert(const_iterator _Where, _Ty&& _Val)
		{	
		return (emplace(_Where, ::std:: move(_Val)));
		}

	template<class... _Valty>
		void emplace_back(_Valty&&... _Val)
		{	
		if (this->_Mylast() == this->_Myend())
			_Reserve(1);
		_Orphan_range(this->_Mylast(), this->_Mylast());
		this->_Getal().construct(_Unfancy(this->_Mylast()),
			::std:: forward<_Valty>(_Val)...);
		++this->_Mylast();
		}

	template<class... _Valty>
		iterator emplace(const_iterator _Where, _Valty&&... _Val)
		{	
		size_type _Off = (_Where)._Ptr - this->_Myfirst();

 




		emplace_back(::std:: forward<_Valty>(_Val)...);
		::std:: rotate(begin() + _Off, end() - 1, end());
		return (begin() + _Off);
		}


	vector(::std:: initializer_list<value_type> _Ilist,
		const _Alloc& _Al = allocator_type())
		: _Mybase(_Al)
		{	
		_Construct(_Ilist.begin(), _Ilist.end());
		}

	_Myt& operator=(::std:: initializer_list<value_type> _Ilist)
		{	
		assign(_Ilist.begin(), _Ilist.end());
		return (*this);
		}

	void assign(::std:: initializer_list<value_type> _Ilist)
		{	
		assign(_Ilist.begin(), _Ilist.end());
		}

	iterator insert(const_iterator _Where,
		::std:: initializer_list<value_type> _Ilist)
		{	
		return (insert(_Where, _Ilist.begin(), _Ilist.end()));
		}

	~vector() noexcept
		{	
		_Tidy();
		}

	_Myt& operator=(const _Myt& _Right)
		{	
		if (this != &_Right)
			{	
			if (this->_Getal() != _Right._Getal()
				&& _Alty::propagate_on_container_copy_assignment::value)
				{	
				_Tidy();
				this->_Copy_alloc(_Right._Getal());
				}

			this->_Orphan_all();

			if (_Right.empty())
				clear();	
			else if (_Right.size() <= size())
				{	
				pointer _Ptr = _Copy_unchecked(_Right._Myfirst(),
					_Right._Mylast(), this->_Myfirst());	
				_Destroy(_Ptr, this->_Mylast());	
				this->_Mylast() = this->_Myfirst() + _Right.size();
				}
			else if (_Right.size() <= capacity())
				{	
				pointer _Ptr = _Right._Myfirst() + size();
				_Copy_unchecked(_Right._Myfirst(),
					_Ptr, this->_Myfirst());
				this->_Mylast() = _Ucopy(_Ptr, _Right._Mylast(),
					this->_Mylast());
				}
			else
				{	
				if (this->_Myfirst() != pointer())
					{	
					_Destroy(this->_Myfirst(), this->_Mylast());
					this->_Getal().deallocate(this->_Myfirst(),
						this->_Myend() - this->_Myfirst());
					}
				if (_Buy(_Right.size()))
					try {
					this->_Mylast() =
						_Ucopy(_Right._Myfirst(), _Right._Mylast(),
						this->_Myfirst());
					} catch (...) {
					_Tidy();
					throw;
					}
				}
			}
		return (*this);
		}

	void reserve(size_type _Count)
		{	
		if (capacity() < _Count)
			{	
			if (max_size() < _Count)
				_Xlen();
			_Reallocate(_Count);
			}
		}

	size_type capacity() const noexcept
		{	
		return (this->_Myend() - this->_Myfirst());
		}

	size_type _Unused_capacity() const noexcept
		{	
		return (this->_Myend() - this->_Mylast());
		}

	size_type _Has_unused_capacity() const noexcept
		{	
		return (this->_Myend() != this->_Mylast());
		}

	iterator begin() noexcept
		{	
		return (iterator(this->_Myfirst(), &this->_Get_data()));
		}

	const_iterator begin() const noexcept
		{	
		return (const_iterator(this->_Myfirst(), &this->_Get_data()));
		}

	iterator end() noexcept
		{	
		return (iterator(this->_Mylast(), &this->_Get_data()));
		}

	const_iterator end() const noexcept
		{	
		return (const_iterator(this->_Mylast(), &this->_Get_data()));
		}

	iterator _Make_iter(const_iterator _Where) const
		{	
		return (iterator(_Where._Ptr, &this->_Get_data()));
		}

	reverse_iterator rbegin() noexcept
		{	
		return (reverse_iterator(end()));
		}

	const_reverse_iterator rbegin() const noexcept
		{	
		return (const_reverse_iterator(end()));
		}

	reverse_iterator rend() noexcept
		{	
		return (reverse_iterator(begin()));
		}

	const_reverse_iterator rend() const noexcept
		{	
		return (const_reverse_iterator(begin()));
		}

	const_iterator cbegin() const noexcept
		{	
		return (begin());
		}

	const_iterator cend() const noexcept
		{	
		return (end());
		}

	const_reverse_iterator crbegin() const noexcept
		{	
		return (rbegin());
		}

	const_reverse_iterator crend() const noexcept
		{	
		return (rend());
		}

	void shrink_to_fit()
		{	
		if (_Has_unused_capacity())
			{	
			if (empty())
				_Tidy();
			else
				_Reallocate(size());
			}
		}

	void resize(size_type _Newsize)
		{	
		if (_Newsize < size())
			_Pop_back_n(size() - _Newsize);
		else if (size() < _Newsize)
			{	
			_Reserve(_Newsize - size());
			try {
			_Uninitialized_default_fill_n(this->_Mylast(), _Newsize - size(),
				this->_Getal());
			} catch (...) {
			_Tidy();
			throw;
			}
			this->_Mylast() += _Newsize - size();
			}
		}

	void resize(size_type _Newsize, const value_type& _Val)
		{	
		if (_Newsize < size())
			_Pop_back_n(size() - _Newsize);
		else if (size() < _Newsize)
			{	
			const value_type *_Ptr = ::std:: addressof(_Val);

			if (_Inside(_Ptr))
				{	
				const difference_type _Idx = _Ptr
					- _Unfancy(this->_Myfirst());
				_Reserve(_Newsize - size());
				_Ptr = _Unfancy(this->_Myfirst()) + _Idx;
				}
			else
				_Reserve(_Newsize - size());

			try {
			_Ufill(this->_Mylast(), _Newsize - size(), _Ptr);
			} catch (...) {
			_Tidy();
			throw;
			}
			this->_Mylast() += _Newsize - size();
			}
		}

	size_type size() const noexcept
		{	
		return (this->_Mylast() - this->_Myfirst());
		}

	size_type max_size() const noexcept
		{	
		return (this->_Getal().max_size());
		}

	bool empty() const noexcept
		{	
		return (this->_Myfirst() == this->_Mylast());
		}

	_Alloc get_allocator() const noexcept
		{	
		_Alloc _Ret(this->_Getal());
		return (_Ret);
		}

	const_reference at(size_type _Pos) const
		{	
		if (size() <= _Pos)
			_Xran();
		return (*(this->_Myfirst() + _Pos));
		}

	reference at(size_type _Pos)
		{	
		if (size() <= _Pos)
			_Xran();
		return (*(this->_Myfirst() + _Pos));
		}

	const_reference operator[](size_type _Pos) const
		{	
 










		return (*(this->_Myfirst() + _Pos));
		}

	reference operator[](size_type _Pos)
		{	
 










		return (*(this->_Myfirst() + _Pos));
		}

	_Ty * data() noexcept
		{	
		return (_Unfancy(this->_Myfirst()));
		}

	const _Ty * data() const noexcept
		{	
		return (_Unfancy(this->_Myfirst()));
		}

	reference front()
		{	
		return (*begin());
		}

	const_reference front() const
		{	
		return (*begin());
		}

	reference back()
		{	
		return (*(end() - 1));
		}

	const_reference back() const
		{	
		return (*(end() - 1));
		}

	void push_back(const value_type& _Val)
		{	
		if (_Inside(::std:: addressof(_Val)))
			{	
			size_type _Idx = ::std:: addressof(_Val) - _Unfancy(this->_Myfirst());
			if (this->_Mylast() == this->_Myend())
				_Reserve(1);
			_Orphan_range(this->_Mylast(), this->_Mylast());
			this->_Getal().construct(_Unfancy(this->_Mylast()),
				this->_Myfirst()[_Idx]);
			++this->_Mylast();
			}
		else
			{	
			if (this->_Mylast() == this->_Myend())
				_Reserve(1);
			_Orphan_range(this->_Mylast(), this->_Mylast());
			this->_Getal().construct(_Unfancy(this->_Mylast()),
				_Val);
			++this->_Mylast();
			}
		}

 













	void pop_back()
		{	
		this->_Getal().destroy(_Unfancy(this->_Mylast() - 1));
		--this->_Mylast();
		}
 

	template<class _Iter>
		typename enable_if<_Is_iterator<_Iter>::value,
			void>::type
		assign(_Iter _First, _Iter _Last)
		{	
		clear();
		_Assign(_First, _Last, _Iter_cat_t<_Iter>());
		}

	template<class _Iter>
		void _Assign(_Iter _First, _Iter _Last,
			input_iterator_tag)
		{	
		for (; _First != _Last; ++_First)
			emplace_back(*_First);
		}

	template<class _Iter>
		void _Assign(_Iter _First, _Iter _Last,
			forward_iterator_tag)
		{	
		size_type _Newsize = ::std:: distance(_First, _Last);

		if (capacity() < _Newsize)
			{	
			size_type _Newcapacity = _Grow_to(_Newsize);
			_Tidy();
			_Buy(_Newcapacity);
			}

		this->_Mylast() = _Ucopy(_First, _Last, this->_Myfirst());
		}

	void assign(size_type _Count, const value_type& _Val)
		{	
		clear();
		insert(begin(), _Count, _Val);
		}

	iterator insert(const_iterator _Where, const _Ty& _Val)
		{	
		return (_Insert_n(_Where, (size_type)1, _Val));
		}

	iterator insert(const_iterator _Where, size_type _Count,
		const _Ty& _Val)
		{	
		return (_Insert_n(_Where, _Count, _Val));
		}

	template<class _Iter>
		typename enable_if<_Is_iterator<_Iter>::value,
			iterator>::type
		insert(const_iterator _Where, _Iter _First, _Iter _Last)
		{	
		size_type _Off = (_Where)._Ptr - this->_Myfirst();
		_Insert(_Where, _First, _Last, _Iter_cat_t<_Iter>());
		return (begin() + _Off);
		}

	template<class _Iter>
		void _Insert(const_iterator _Where,
			_Iter _First, _Iter _Last,
				input_iterator_tag)
		{	
		size_type _Off = (_Where)._Ptr - this->_Myfirst();

 




		if (_First != _Last)
			{	
			size_type _Oldsize = size();

			try {
			for (; _First != _Last; ++_First)
				push_back(*_First);	

			} catch (...) {
			erase(begin() + _Oldsize, end());
			throw;
			}

			::std:: rotate(begin() + _Off, begin() + _Oldsize, end());
			}
		}

	template<class _Iter>
		void _Insert(const_iterator _Where,
			_Iter _First, _Iter _Last,
				forward_iterator_tag)
		{	
 







		size_type _Count = ::std:: distance(_First, _Last);
		if (_Count == 0)
			;
		else if (_Unused_capacity() < _Count)
			{	
			if (max_size() - size() < _Count)
				_Xlen();	

			size_type _Capacity = _Grow_to(size() + _Count);
			pointer _Newvec = this->_Getal().allocate(_Capacity);
			pointer _Ptr = _Newvec;

			try {
			_Ptr = _Umove(this->_Myfirst(), (_Where)._Ptr,
				_Newvec);	
			_Ptr = _Ucopy(_First, _Last, _Ptr);	
			_Umove((_Where)._Ptr, this->_Mylast(),
				_Ptr);	
			} catch (...) {
			_Destroy(_Newvec, _Ptr);
			this->_Getal().deallocate(_Newvec, _Capacity);
			throw;
			}

			_Count += size();
			if (this->_Myfirst() != pointer())
				{	
				_Destroy(this->_Myfirst(), this->_Mylast());
				this->_Getal().deallocate(this->_Myfirst(),
					this->_Myend() - this->_Myfirst());
				}

			this->_Orphan_all();
			this->_Myend() = _Newvec + _Capacity;
			this->_Mylast() = _Newvec + _Count;
			this->_Myfirst() = _Newvec;
			}
		else
			{	
			_Ucopy(_First, _Last, this->_Mylast());
			::std:: rotate((_Where)._Ptr, this->_Mylast(),
				this->_Mylast() + _Count);
			this->_Mylast() += _Count;
			_Orphan_range((_Where)._Ptr, this->_Mylast());
			}
		}

 














	iterator erase(const_iterator _Where)
		{	
		_Move_unchecked((_Where)._Ptr + 1, this->_Mylast(),
			(_Where)._Ptr);
		_Destroy(this->_Mylast() - 1, this->_Mylast());
		--this->_Mylast();
		return (_Make_iter(_Where));
		}
 

	iterator erase(const_iterator _First_arg,
		const_iterator _Last_arg)
		{	
		if (_First_arg == begin() && _Last_arg == end())
			clear();
		else if (_First_arg != _Last_arg)
			{	
			iterator _First = _Make_iter(_First_arg);
			iterator _Last = _Make_iter(_Last_arg);

			if (_First != _Last)
				{	
 









				pointer _Ptr = _Move_unchecked((_Last)._Ptr, this->_Mylast(),
					(_First)._Ptr);
 

				_Destroy(_Ptr, this->_Mylast());
				this->_Mylast() = _Ptr;
				}
			}
		return (_Make_iter(_First_arg));
		}

	void _Pop_back_n(size_type _Count)
		{	
		pointer _Ptr = this->_Mylast() - _Count;

 



		_Destroy(_Ptr, this->_Mylast());
		this->_Mylast() = _Ptr;
		}

	void clear() noexcept
		{	
		this->_Orphan_all();
		_Destroy(this->_Myfirst(), this->_Mylast());
		this->_Mylast() = this->_Myfirst();
		}

	void swap(_Myt& _Right)
		noexcept(_Alty::propagate_on_container_swap::value || _Alty::is_always_equal::value)
		{	
		if (this != &_Right)
			{	
			_Pocs(this->_Getal(), _Right._Getal());
			this->_Swap_all(_Right);
			_Swap_adl(this->_Myfirst(), _Right._Myfirst());
			_Swap_adl(this->_Mylast(), _Right._Mylast());
			_Swap_adl(this->_Myend(), _Right._Myend());
			}
		}

protected:
	bool _Buy(size_type _Capacity)
		{	
		this->_Myfirst() = pointer();
		this->_Mylast() = pointer();
		this->_Myend() = pointer();

		if (_Capacity == 0)
			return (false);
		else if (max_size() < _Capacity)
			_Xlen();	
		else
			{	
			this->_Myfirst() = this->_Getal().allocate(_Capacity);
			this->_Mylast() = this->_Myfirst();
			this->_Myend() = this->_Myfirst() + _Capacity;
			}
		return (true);
		}

	void _Destroy(pointer _First, pointer _Last)
		{	
		_Destroy_range(_First, _Last, this->_Getal());
		}

	size_type _Grow_to(size_type _Count) const
		{	
		size_type _Capacity = capacity();

		_Capacity = max_size() - _Capacity / 2 < _Capacity
			? 0 : _Capacity + _Capacity / 2;	
		if (_Capacity < _Count)
			_Capacity = _Count;
		return (_Capacity);
		}

	bool _Inside(const value_type *_Ptr) const
		{	
		return (_Ptr < _Unfancy(this->_Mylast()) && _Unfancy(this->_Myfirst()) <= _Ptr);
		}

	void _Reallocate(size_type _Count)
		{	
		pointer _Ptr = this->_Getal().allocate(_Count);

		try {
		_Umove(this->_Myfirst(), this->_Mylast(), _Ptr);
		} catch (...) {
		this->_Getal().deallocate(_Ptr, _Count);
		throw;
		}

		size_type _Size = size();
		if (this->_Myfirst() != pointer())
			{	
			_Destroy(this->_Myfirst(), this->_Mylast());
			this->_Getal().deallocate(this->_Myfirst(),
				this->_Myend() - this->_Myfirst());
			}

		this->_Orphan_all();
		this->_Myend() = _Ptr + _Count;
		this->_Mylast() = _Ptr + _Size;
		this->_Myfirst() = _Ptr;
		}

	void _Reserve(size_type _Count)
		{	
		if (_Unused_capacity() < _Count)
			{	
			if (max_size() - size() < _Count)
				_Xlen();
			_Reallocate(_Grow_to(size() + _Count));
			}
		}

	void _Tidy()
		{	
		if (this->_Myfirst() != pointer())
			{	
			this->_Orphan_all();
			_Destroy(this->_Myfirst(), this->_Mylast());
			this->_Getal().deallocate(this->_Myfirst(),
				this->_Myend() - this->_Myfirst());
			this->_Myfirst() = pointer();
			this->_Mylast() = pointer();
			this->_Myend() = pointer();
			}
		}

	template<class _Iter>
		pointer _Ucopy(_Iter _First, _Iter _Last, pointer _Ptr)
		{	
		return (_Uninitialized_copy(_First, _Last,
			_Ptr, this->_Getal()));
		}

	template<class _Iter>
		pointer _Umove(_Iter _First, _Iter _Last, pointer _Ptr)
		{	
		return (_Uninitialized_move(_First, _Last,
			_Ptr, this->_Getal()));
		}

	iterator _Insert_n(const_iterator _Where,
		size_type _Count, const value_type& _Val)
		{	
 






		size_type _Off = (_Where)._Ptr - this->_Myfirst();
		if (_Count == 0)
			;
		else if (_Unused_capacity() < _Count)
			{	
			if (max_size() - size() < _Count)
				_Xlen();	

			size_type _Capacity = _Grow_to(size() + _Count);
			pointer _Newvec = this->_Getal().allocate(_Capacity);
			size_type _Whereoff = (_Where)._Ptr - this->_Myfirst();
			int _Ncopied = 0;

			try {
			_Ufill(_Newvec + _Whereoff, _Count,
				::std:: addressof(_Val));	
			++_Ncopied;
			_Umove(this->_Myfirst(), (_Where)._Ptr,
				_Newvec);	
			++_Ncopied;
			_Umove((_Where)._Ptr, this->_Mylast(),
				_Newvec + (_Whereoff + _Count));	
			} catch (...) {
			if (1 < _Ncopied)
				_Destroy(_Newvec, _Newvec + _Whereoff);
			if (0 < _Ncopied)
				_Destroy(_Newvec + _Whereoff, _Newvec + _Whereoff + _Count);
			this->_Getal().deallocate(_Newvec, _Capacity);
			throw;
			}

			_Count += size();
			if (this->_Myfirst() != pointer())
				{	
				_Destroy(this->_Myfirst(), this->_Mylast());
				this->_Getal().deallocate(this->_Myfirst(),
					this->_Myend() - this->_Myfirst());
				}

			this->_Orphan_all();
			this->_Myend() = _Newvec + _Capacity;
			this->_Mylast() = _Newvec + _Count;
			this->_Myfirst() = _Newvec;
			}
		else if ((size_type)(this->_Mylast() - (_Where)._Ptr)
			< _Count)
			{	
			value_type _Tmp = _Val;	

			_Umove((_Where)._Ptr, this->_Mylast(),
				(_Where)._Ptr + _Count);	

			try {
			_Ufill(this->_Mylast(),
				_Count - (this->_Mylast() - (_Where)._Ptr),
				::std:: addressof(_Tmp));	
			} catch (...) {
			_Destroy((_Where)._Ptr + _Count,
				this->_Mylast() + _Count);
			throw;
			}

			this->_Mylast() += _Count;
			_Orphan_range((_Where)._Ptr, this->_Mylast());
			::std:: fill((_Where)._Ptr, this->_Mylast() - _Count,
				_Tmp);	
			}
		else
			{	
			value_type _Tmp = _Val;	

			pointer _Oldend = this->_Mylast();
			this->_Mylast() = _Umove(_Oldend - _Count, _Oldend,
				this->_Mylast());	

			_Orphan_range((_Where)._Ptr, this->_Mylast());
			_Move_backward_unchecked((_Where)._Ptr, _Oldend - _Count,
				_Oldend);	
			::std:: fill((_Where)._Ptr,
				(_Where)._Ptr + _Count, _Tmp);	
			}
		return (begin() + _Off);
		}

	pointer _Ufill(pointer _Ptr, size_type _Count, const value_type *_Pval)
		{	
		_Uninitialized_fill_n(_Ptr, _Count, _Pval, this->_Getal());
		return (_Ptr + _Count);
		}

	[[noreturn]] void _Xlen() const
		{	
		_Xlength_error("vector<T> too long");
		}

	[[noreturn]] void _Xran() const
		{	
		_Xout_of_range("invalid vector<T> subscript");
		}

 


















	void _Orphan_range(pointer, pointer) const
		{	
		}
 
	};

		

template<class _Ty,
	class _Alloc> inline
	void swap(vector<_Ty, _Alloc>& _Left, vector<_Ty, _Alloc>& _Right)
		noexcept(noexcept(_Left.swap(_Right)))
	{	
	_Left.swap(_Right);
	}

template<class _Ty,
	class _Alloc> inline
	bool operator==(const vector<_Ty, _Alloc>& _Left,
		const vector<_Ty, _Alloc>& _Right)
	{	
	return (_Left.size() == _Right.size()
		&& ::std:: equal(_Left.begin(), _Left.end(), _Right.begin()));
	}

template<class _Ty,
	class _Alloc> inline
	bool operator!=(const vector<_Ty, _Alloc>& _Left,
		const vector<_Ty, _Alloc>& _Right)
	{	
	return (!(_Left == _Right));
	}

template<class _Ty,
	class _Alloc> inline
	bool operator<(const vector<_Ty, _Alloc>& _Left,
		const vector<_Ty, _Alloc>& _Right)
	{	
	return (::std:: lexicographical_compare(_Left.begin(), _Left.end(),
		_Right.begin(), _Right.end()));
	}

template<class _Ty,
	class _Alloc> inline
	bool operator>(const vector<_Ty, _Alloc>& _Left,
		const vector<_Ty, _Alloc>& _Right)
	{	
	return (_Right < _Left);
	}

template<class _Ty,
	class _Alloc> inline
	bool operator<=(const vector<_Ty, _Alloc>& _Left,
		const vector<_Ty, _Alloc>& _Right)
	{	
	return (!(_Right < _Left));
	}

template<class _Ty,
	class _Alloc> inline
	bool operator>=(const vector<_Ty, _Alloc>& _Left,
		const vector<_Ty, _Alloc>& _Right)
	{	
	return (!(_Left < _Right));
	}




typedef unsigned int _Vbase;	
const int _VBITS = 8 * sizeof (_Vbase);	

		
template<class _Alloc>
	class _Vb_iter_base
		: public _Iterator012<random_access_iterator_tag,
			bool,
			typename _Alloc::difference_type,
			bool *,
			bool,
			_Iterator_base>
	{	
public:
	typedef typename _Alloc::size_type _Sizet;
	typedef vector<bool, _Alloc> _Mycont;

	_Vb_iter_base()
		: _Myptr(0), _Myoff(0)
		{	
		}

	_Vb_iter_base(const _Vbase *_Ptr, _Sizet _Off,
		const _Container_base *_Mypvbool)
		: _Myptr(_Ptr), _Myoff(_Off)
		{	
		this->_Adopt(_Mypvbool);
		}

	void _Advance(_Sizet _Off)
		{	
		_Myoff += _Off;
		_Myptr += _Myoff / _VBITS;
		_Myoff %= _VBITS;
		}

	int _Valid(_Sizet _Inc) const
		{	
 








		(void) _Inc;
		return (-1);
 
		}

	const _Vbase *_Myptr;
	_Sizet _Myoff;
	};

		
template<class _Alloc>
	class _Vb_reference
		: public _Vb_iter_base<_Alloc>
	{	
	typedef _Vb_iter_base<_Alloc> _Mybase;
	typedef _Vb_reference<_Alloc> _Mytype;

	_Vb_reference() noexcept
		{	
		}

public:
	_Vb_reference(const _Mybase& _Right)
		: _Mybase(_Right._Myptr, _Right._Myoff, _Right._Getcont())
		{	
		}

	_Mytype& operator=(const _Mytype& _Right) noexcept
		{	
		return (*this = bool(_Right));
		}

	_Mytype& operator=(bool _Val) noexcept
		{	
		if (_Val)
			*(_Vbase *)_Getptr() |= _Mask();
		else
			*(_Vbase *)_Getptr() &= (~_Mask());	
		return (*this);
		}

	void flip() noexcept
		{	
		*(_Vbase *)_Getptr() ^= _Mask();
		}

	operator bool() const noexcept
		{	
		return ((*_Getptr() & _Mask()) != 0);
		}

	const _Vbase *_Getptr() const
		{	
 













		return (this->_Myptr);
		}

protected:
	_Vbase _Mask() const
		{	
		return ((_Vbase)(1) << this->_Myoff);
		}
	};

template<class _Alloc> inline
	void swap(_Vb_reference<_Alloc> _Left,
		_Vb_reference<_Alloc> _Right)
	{	
	bool _Val = _Left;	
	_Left = _Right;
	_Right = _Val;
	}

		
template<class _Alloc>
	class _Vb_const_iterator
		: public _Vb_iter_base<_Alloc>
	{	
public:
	typedef _Vb_iter_base<_Alloc> _Mybase;
	typedef _Vb_const_iterator<_Alloc> _Mytype;

	typedef _Vb_reference<_Alloc> _Reft;
	typedef bool const_reference;

	typedef random_access_iterator_tag iterator_category;
	typedef bool value_type;
	typedef typename _Alloc::size_type size_type;
	typedef typename _Alloc::difference_type difference_type;
	typedef const_reference *pointer;
	typedef const_reference reference;

	_Vb_const_iterator()
		{	
		}

	_Vb_const_iterator(const _Vbase *_Ptr, const _Container_base *_Mypvbool)
		: _Mybase(_Ptr, 0, _Mypvbool)
		{	
		}

	const_reference operator*() const
		{	
		return (_Reft(*this));
		}

	_Mytype& operator++()
		{	
		_Inc();
		return (*this);
		}

	_Mytype operator++(int)
		{	
		_Mytype _Tmp = *this;
		++*this;
		return (_Tmp);
		}

	_Mytype& operator--()
		{	
		_Dec();
		return (*this);
		}

	_Mytype operator--(int)
		{	
		_Mytype _Tmp = *this;
		--*this;
		return (_Tmp);
		}

	_Mytype& operator+=(difference_type _Off)
		{	
		if (_Off < 0 && this->_Myoff < 0 - (size_type)_Off)
			{	
			this->_Myoff += _Off;
			this->_Myptr -= 1 + ((size_type)(-1) - this->_Myoff) / _VBITS;
			this->_Myoff %= _VBITS;
			}
		else
			{	
			this->_Myoff += _Off;
			this->_Myptr += this->_Myoff / _VBITS;
			this->_Myoff %= _VBITS;
			}
		return (*this);
		}

	_Mytype operator+(difference_type _Off) const
		{	
		_Mytype _Tmp = *this;
		return (_Tmp += _Off);
		}

	_Mytype& operator-=(difference_type _Off)
		{	
		return (*this += -_Off);
		}

	_Mytype operator-(difference_type _Off) const
		{	
		_Mytype _Tmp = *this;
		return (_Tmp -= _Off);
		}

	difference_type operator-(
		const _Mytype& _Right) const
		{	
		_Compat(_Right);
		return (_VBITS * (this->_Myptr - _Right._Myptr)
			+ (difference_type)this->_Myoff
			- (difference_type)_Right._Myoff);
		}

	const_reference operator[](difference_type _Off) const
		{	
		return (*(*this + _Off));
		}

	bool operator==(const _Mytype& _Right) const
		{	
		_Compat(_Right);
		return (this->_Myptr == _Right._Myptr
			&& this->_Myoff == _Right._Myoff);
		}

	bool operator!=(const _Mytype& _Right) const
		{	
		return (!(*this == _Right));
		}

	bool operator<(const _Mytype& _Right) const
		{	
		_Compat(_Right);
		return (this->_Myptr < _Right._Myptr
			|| (this->_Myptr == _Right._Myptr
				&& this->_Myoff < _Right._Myoff));
		}

	bool operator>(const _Mytype& _Right) const
		{	
		return (_Right < *this);
		}

	bool operator<=(const _Mytype& _Right) const
		{	
		return (!(_Right < *this));
		}

	bool operator>=(const _Mytype& _Right) const
		{	
		return (!(*this < _Right));
		}

 













	void _Compat(const _Mytype&) const
		{	
		}
 

	void _Dec()
		{	
		if (this->_Myoff != 0)
			--this->_Myoff;
		else
			{	
 











			this->_Myoff = _VBITS - 1;
			--this->_Myptr;
			}
		}

	void _Inc()
		{	
		if (this->_Myoff < _VBITS - 1)
			++this->_Myoff;
		else
			{	
 











			this->_Myoff = 0;
			++this->_Myptr;
			}
		}
	};

template<class _Alloc> inline
	_Vb_const_iterator<_Alloc> operator+(
		typename _Alloc::difference_type _Off,
		_Vb_const_iterator<_Alloc> _Right)
		{	
		return (_Right += _Off);
		}

template<class _Alloc>
	struct _Is_checked_helper<_Vb_const_iterator<_Alloc> >
		: public true_type
	{	
	};

	
template<class _Alloc>
	class _Vb_iterator
		: public _Vb_const_iterator<_Alloc>
	{	
public:
	typedef _Vb_const_iterator<_Alloc> _Mybase;
	typedef _Vb_iterator<_Alloc> _Mytype;

	typedef _Vb_reference<_Alloc> _Reft;
	typedef bool const_reference;

	typedef random_access_iterator_tag iterator_category;
	typedef bool value_type;
	typedef typename _Alloc::size_type size_type;
	typedef typename _Alloc::difference_type difference_type;
	typedef _Reft *pointer;
	typedef _Reft reference;

	_Vb_iterator()
		{	
		}

	_Vb_iterator(_Vbase *_Ptr, _Container_base *_Mypvbool)
		: _Mybase(_Ptr, _Mypvbool)
		{	
		}

	reference operator*() const
		{	
		return (_Reft(*this));
		}

	_Mytype& operator++()
		{	
		++*(_Mybase *)this;
		return (*this);
		}

	_Mytype operator++(int)
		{	
		_Mytype _Tmp = *this;
		++*this;
		return (_Tmp);
		}

	_Mytype& operator--()
		{	
		--*(_Mybase *)this;
		return (*this);
		}

	_Mytype operator--(int)
		{	
		_Mytype _Tmp = *this;
		--*this;
		return (_Tmp);
		}

	_Mytype& operator+=(difference_type _Off)
		{	
		*(_Mybase *)this += _Off;
		return (*this);
		}

	_Mytype operator+(difference_type _Off) const
		{	
		_Mytype _Tmp = *this;
		return (_Tmp += _Off);
		}

	_Mytype& operator-=(difference_type _Off)
		{	
		return (*this += -_Off);
		}

	_Mytype operator-(difference_type _Off) const
		{	
		_Mytype _Tmp = *this;
		return (_Tmp -= _Off);
		}

	difference_type operator-(const _Mybase& _Right) const
		{	
		return (*(_Mybase *)this - _Right);
		}

	reference operator[](difference_type _Off) const
		{	
		return (*(*this + _Off));
		}
	};

template<class _Alloc> inline
	_Vb_iterator<_Alloc> operator+(typename _Alloc::difference_type _Off,
		_Vb_iterator<_Alloc> _Right)
		{	
		return (_Right += _Off);
		}

template<class _Alloc>
	struct _Is_checked_helper<_Vb_iterator<_Alloc> >
		: public true_type
	{	
	};

		
template<class _Alloc>
	class _Vb_val
		: public _Container_base
	{	
public:
	typedef vector<_Vbase, _Alloc> _Vectype;
	typedef typename _Vectype::_Alty _Alty;
	typedef typename _Alty::size_type size_type;

	_Vb_val(size_type _Count, const bool& _Val)
		: _Myvec(_Nw(_Count), (_Vbase) (_Val ? -1 : 0))
		{	
		_Alloc_proxy();
		_Mysize = 0;
		}

	_Vb_val(size_type _Count, const bool& _Val, const _Alloc& _Al)
		: _Myvec(_Nw(_Count), (_Vbase)(_Val ? -1 : 0), _Al)
		{	
		_Alloc_proxy();
		_Mysize = 0;
		}

	_Vb_val(const _Vb_val& _Right)
		: _Myvec(_Right._Myvec),
			_Mysize(_Right._Mysize)
		{	
		_Alloc_proxy();
		}

	_Vb_val(const _Vb_val& _Right, const _Alloc& _Al)
		: _Myvec(_Right._Myvec, _Al),
			_Mysize(_Right._Mysize)
		{	
		_Alloc_proxy();
		}

	_Vb_val(_Vb_val&& _Right)
		: _Myvec(::std:: forward<_Vectype>(_Right._Myvec)),
			_Mysize(_Right._Mysize)
		{	
		_Right._Mysize = 0;
		_Alloc_proxy();
		}

	_Vb_val(_Vb_val&& _Right, const _Alloc& _Al)
		: _Myvec(::std:: forward<_Vectype>(_Right._Myvec), _Al),
			_Mysize(_Right._Mysize)
		{	
		_Right._Mysize = 0;
		_Alloc_proxy();
		}

	~_Vb_val() noexcept
		{	
		_Free_proxy();
		}

 
	void _Alloc_proxy()
		{	
		}

	void _Free_proxy()
		{	
		}

 




















	static size_type _Nw(size_type _Count)
		{	
		return ((_Count + _VBITS - 1) / _VBITS);
		}

	_Vectype _Myvec;	
	typename _Alty::size_type _Mysize;	
	};

		

template<class _Alloc>
	class vector<bool, _Alloc>
		: public _Vb_val<_Alloc>
	{	
public:
	typedef vector<bool, _Alloc> _Myt;
	typedef _Vb_val<_Alloc> _Mybase;
	typedef typename _Mybase::_Alty _Alty;
	typedef typename _Mybase::_Vectype _Vectype;

	typedef typename _Alty::size_type size_type;
	typedef typename _Alty::difference_type difference_type;
	typedef bool _Ty;
	typedef _Alloc allocator_type;

	typedef _Vb_reference<_Alty> reference;
	typedef bool const_reference;
	typedef bool value_type;

	typedef reference _Reft;
	typedef _Vb_const_iterator<_Alty> const_iterator;
	typedef _Vb_iterator<_Alty> iterator;

	typedef iterator pointer;
	typedef const_iterator const_pointer;
	typedef ::std:: reverse_iterator<iterator> reverse_iterator;
	typedef ::std:: reverse_iterator<const_iterator> const_reverse_iterator;

	static const int _VBITS = ::std:: _VBITS;
	enum {_EEN_VBITS = _VBITS};	
	vector()
		: _Mybase(0, false)
		{	
		}

	explicit vector(const _Alloc& _Al)
		: _Mybase(0, false, _Al)
		{	
		}

	explicit vector(size_type _Count, const _Alloc& _Al = _Alloc())
		: _Mybase(_Count, false, _Al)
		{	
		_Trim(_Count);
		}

	vector(size_type _Count, const bool& _Val, const _Alloc& _Al = _Alloc())
		: _Mybase(_Count, _Val, _Al)
		{	
		_Trim(_Count);
		}

	vector(const _Myt& _Right)
		: _Mybase(_Right)
		{	
		}

	vector(const _Myt& _Right, const _Alloc& _Al)
		: _Mybase(_Right, _Al)
		{	
		}

	template<class _Iter,
		class = typename enable_if<_Is_iterator<_Iter>::value,
			void>::type>
		vector(_Iter _First, _Iter _Last, const _Alloc& _Al = _Alloc())
		: _Mybase(0, false, _Al)
		{	
		_BConstruct(_First, _Last);
		}

	template<class _Iter>
		void _BConstruct(_Iter _First, _Iter _Last)
		{	
		insert(begin(), _First, _Last);
		}

	vector(_Myt&& _Right)
		: _Mybase(::std:: forward<_Myt>(_Right))
		{	
		}

	vector(_Myt&& _Right, const _Alloc& _Al)
		: _Mybase(::std:: forward<_Myt>(_Right), _Al)
		{	
		}

	_Myt& operator=(_Myt&& _Right)
		{	
		if (this != &_Right)
			{	
			clear();

			if (_Alty::propagate_on_container_move_assignment::value
				&& this->get_allocator() != _Right.get_allocator())
				{	
				this->_Free_proxy();
				this->_Myvec = ::std:: move(_Right._Myvec);
				this->_Alloc_proxy();
				}
			else
				this->_Myvec = ::std:: move(_Right._Myvec);


			this->_Mysize = _Right._Mysize;
			_Right._Mysize = 0;
			}
		return (*this);
		}

	template<class... _Valty>
		void emplace_back(_Valty&&... _Val)
		{	
		bool _Tmp(::std:: forward<_Valty>(_Val)...);
		push_back(_Tmp);
		}

	template<class... _Valty>
		iterator emplace(const_iterator _Where, _Valty&&... _Val)
		{	
		bool _Tmp(::std:: forward<_Valty>(_Val)...);
		return (insert(_Where, _Tmp));
		}


	vector(::std:: initializer_list<bool> _Ilist,
			const _Alloc& _Al = allocator_type())
		: _Mybase(0, false, _Al)
		{	
		insert(begin(), _Ilist.begin(), _Ilist.end());
		}

	_Myt& operator=(::std:: initializer_list<bool> _Ilist)
		{	
		assign(_Ilist.begin(), _Ilist.end());
		return (*this);
		}

	void assign(::std:: initializer_list<bool> _Ilist)
		{	
		assign(_Ilist.begin(), _Ilist.end());
		}

	iterator insert(const_iterator _Where,
			::std:: initializer_list<bool> _Ilist)
		{	
		return (insert(_Where, _Ilist.begin(), _Ilist.end()));
		}

	~vector() noexcept
		{	
		}

	_Myt& operator=(const _Myt& _Right)
		{	
		this->_Mysize = _Right._Mysize;
		this->_Myvec = _Right._Myvec;
		return (*this);
		}

	void reserve(size_type _Count)
		{	
		this->_Myvec.reserve(this->_Nw(_Count));
		}

	size_type capacity() const noexcept
		{	
		return (this->_Myvec.capacity() * _VBITS);
		}

	iterator begin() noexcept
		{	
		return (iterator(this->_Myvec.data(), this));
		}

	const_iterator begin() const noexcept
		{	
		return (const_iterator(this->_Myvec.data(), this));
		}

	iterator end() noexcept
		{	
		iterator _Tmp = begin();
		if (0 < this->_Mysize)
			_Tmp += this->_Mysize;
		return (_Tmp);
		}

	const_iterator end() const noexcept
		{	
		const_iterator _Tmp = begin();
		if (0 < this->_Mysize)
			_Tmp += this->_Mysize;
		return (_Tmp);
		}

	const_iterator cbegin() const noexcept
		{	
		return (begin());
		}

	const_iterator cend() const noexcept
		{	
		return (end());
		}

	const_reverse_iterator crbegin() const noexcept
		{	
		return (rbegin());
		}

	const_reverse_iterator crend() const noexcept
		{	
		return (rend());
		}

	void shrink_to_fit()
		{	
		if (this->_Myvec._Has_unused_capacity())
			{	
			_Myt _Tmp(*this);
			swap(_Tmp);
			}
		}

	iterator _Make_iter(const_iterator _Where)
		{	
		iterator _Tmp = begin();
		if (0 < this->_Mysize)
			_Tmp += _Where - begin();
		return (_Tmp);
		}

	reverse_iterator rbegin() noexcept
		{	
		return (reverse_iterator(end()));
		}

	const_reverse_iterator rbegin() const noexcept
		{	
		return (const_reverse_iterator(end()));
		}

	reverse_iterator rend() noexcept
		{	
		return (reverse_iterator(begin()));
		}

	const_reverse_iterator rend() const noexcept
		{	
		return (const_reverse_iterator(begin()));
		}

	void resize(size_type _Newsize, bool _Val = false)
		{	
		if (size() < _Newsize)
			_Insert_n(end(), _Newsize - size(), _Val);
		else if (_Newsize < size())
			erase(begin() + _Newsize, end());
		}

	size_type size() const noexcept
		{	
		return (this->_Mysize);
		}

	size_type max_size() const noexcept
		{	
		const size_type _Maxsize = this->_Myvec.max_size();
		return (_Maxsize < (size_type)(-1) / _VBITS
			? _Maxsize * _VBITS : (size_type)(-1));
		}

	bool empty() const noexcept
		{	
		return (size() == 0);
		}

	_Alloc get_allocator() const noexcept
		{	
		_Alloc _Ret(this->_Myvec.get_allocator());
		return (_Ret);
		}

	const_reference at(size_type _Off) const
		{	
		if (size() <= _Off)
			_Xran();
		return ((*this)[_Off]);
		}

	reference at(size_type _Off)
		{	
		if (size() <= _Off)
			_Xran();
		return ((*this)[_Off]);
		}

	const_reference operator[](size_type _Off) const
		{	
		const_iterator _It = begin();
		_It._Advance(_Off);
		return (*_It);
		}

	reference operator[](size_type _Off)
		{	
		iterator _It = begin();
		_It._Advance(_Off);
		return (*_It);
		}

	reference front()
		{	
		return (*begin());
		}

	const_reference front() const
		{	
		return (*begin());
		}

	reference back()
		{	
		return (*(end() - 1));
		}

	const_reference back() const
		{	
		return (*(end() - 1));
		}

	void push_back(const bool& _Val)
		{	
		insert(end(), _Val);
		}

	void pop_back()
		{	
		erase(end() - 1);
		}

	template<class _Iter>
		typename enable_if<_Is_iterator<_Iter>::value,
			void>::type
		assign(_Iter _First, _Iter _Last)
		{	
		erase(begin(), end());
		insert(begin(), _First, _Last);
		}

	void assign(size_type _Count, const bool& _Val)
		{	
		erase(begin(), end());
		_Insert_n(begin(), _Count, _Val);
		}

	iterator insert(const_iterator _Where, const bool& _Val)
		{	
		return (_Insert_n(_Where, (size_type)1, _Val));
		}

	iterator insert(const_iterator _Where, size_type _Count,
		const bool& _Val)
		{	
		return (_Insert_n(_Where, _Count, _Val));
		}

	template<class _Iter>
		typename enable_if<_Is_iterator<_Iter>::value,
			iterator>::type
		insert(const_iterator _Where, _Iter _First, _Iter _Last)
		{	
		size_type _Off = _Where - begin();
		_Insert(_Where, _First, _Last, _Iter_cat_t<_Iter>());
		return (begin() + _Off);
		}

	template<class _Iter>
		void _Insert(const_iterator _Where,
			_Iter _First, _Iter _Last,
				input_iterator_tag)
		{	
		size_type _Off = _Where - begin();

		for (; _First != _Last; ++_First, (void)++_Off)
			insert(begin() + _Off, *_First);
		}

	template<class _Iter>
		void _Insert(const_iterator _Where,
			_Iter _First, _Iter _Last,
			forward_iterator_tag)
		{	
		;
		size_type _Count = ::std:: distance(_First, _Last);
		size_type _Off = _Insert_x(_Where, _Count);
		::std:: copy(_First, _Last, begin() + _Off);
		}

	iterator erase(const_iterator _Where_arg)
		{	
		iterator _Where = _Make_iter(_Where_arg);
		size_type _Off = _Where - begin();

 






		::std:: copy(_Where + 1, end(), _Where);
 

		_Trim(this->_Mysize - 1);
		return (begin() + _Off);
		}

	iterator erase(const_iterator _First_arg,
		const_iterator _Last_arg)
		{	
		iterator _First = _Make_iter(_First_arg);
		iterator _Last = _Make_iter(_Last_arg);
		size_type _Off = _First - begin();

		if (_First != _Last)
			{	
 








			iterator _Next = ::std:: copy(_Last, end(), _First);
			_Trim(_Next - begin());
 
			}
		return (begin() + _Off);
		}

	void clear() noexcept
		{	
		erase(begin(), end());
		}

	void flip() noexcept
		{	
		for (typename _Vectype::iterator _Next = this->_Myvec.begin();
			_Next != this->_Myvec.end(); ++_Next)
			*_Next = (_Vbase)~*_Next;
		_Trim(this->_Mysize);
		}

	void swap(_Myt& _Right)
		{	
		if (this != &_Right)
			{	
			this->_Swap_all(_Right);
			this->_Myvec.swap(_Right._Myvec);
			::std:: swap(this->_Mysize, _Right._Mysize);
			}
		}

	static void swap(reference _Left, reference _Right) noexcept
		{	
		bool _Val = _Left;	

		_Left = _Right;
		_Right = _Val;
		}

	size_t hash() const
		{	
		return (_Hash_seq((const unsigned char *)this->_Myvec.data(),
			this->_Myvec.size() * sizeof (_Vbase)));
		}

	iterator _Insert_n(const_iterator _Where,
		size_type _Count, const bool& _Val)
		{	
		size_type _Off = _Insert_x(_Where, _Count);
		::std:: fill(begin() + _Off, begin() + (_Off + _Count), _Val);
		return (begin() + _Off);
		}

	size_type _Insert_x(const_iterator _Where, size_type _Count)
		{	
		size_type _Off = _Where - begin();

 





		if (_Count == 0)
			;
		else if (max_size() - size() < _Count)
			_Xlen();	
		else
			{	
			this->_Myvec.resize(this->_Nw(size() + _Count), 0);
			if (empty())
				this->_Mysize += _Count;
			else
				{	
				iterator _Oldend = end();
				this->_Mysize += _Count;
				::std:: copy_backward(begin() + _Off, _Oldend, end());
				}

 


			}
		return (_Off);
		}

 
























	void _Orphan_range(size_type, size_type) const
		{	
		}
 

	void _Trim(size_type _Size)
		{	
		if (max_size() < _Size)
			_Xlen();	
		size_type _Words = this->_Nw(_Size);

		if (_Words < this->_Myvec.size())
			this->_Myvec.erase(this->_Myvec.begin() + _Words,
				this->_Myvec.end());
		this->_Mysize = _Size;
		_Size %= _VBITS;
		if (0 < _Size)
			this->_Myvec[_Words - 1] &= ((_Vbase)(1) << _Size) - 1;
		}

	[[noreturn]] void _Xlen() const
		{	
		_Xlength_error("vector<bool> too long");
		}

	[[noreturn]] void _Xran() const
		{	
		_Xout_of_range("invalid vector<bool> subscript");
		}
	};

template<class _Alloc> inline
	bool operator==(const vector<bool, _Alloc>& _Left,
		const vector<bool, _Alloc>& _Right)
	{	
	return (_Left.size() == _Right.size()
		&& ::std:: equal(_Left._Myvec.begin(), _Left._Myvec.end(),
			_Right._Myvec.begin()));
	}

template<class _Alloc> inline
	bool operator!=(const vector<bool, _Alloc>& _Left,
		const vector<bool, _Alloc>& _Right)
	{	
	return (!(_Left == _Right));
	}

	
template<class _Alloc>
	struct hash<vector<bool, _Alloc> >
	{	
	typedef vector<bool, _Alloc> argument_type;
	typedef size_t result_type;

	size_t operator()(const argument_type& _Keyval) const
		{	
		return (_Keyval.hash());
		}
	};
}

 
 #pragma warning(pop)
 #pragma pack(pop)











#pragma once






#pragma once







 #pragma pack(push,8)
 #pragma warning(push,3)
 
 

namespace std {
	
template<bool _Same,
	class _Dest,
	class... _Srcs>
	struct _Tuple_implicit_val0
		: false_type
	{	
	};

template<class... _Dests,
	class... _Srcs>
	struct _Tuple_implicit_val0<true, tuple<_Dests...>, _Srcs...>
		: conjunction<
			is_constructible<_Dests, _Srcs>...,
			is_convertible<_Srcs, _Dests>...
		>::type
	{	
	};

template<class _Dest,
	class... _Srcs>
	struct _Tuple_implicit_val
		: _Tuple_implicit_val0<tuple_size<_Dest>::value == sizeof...(_Srcs), _Dest, _Srcs...>::type
	{	
	};

template<class _Dest,
	class... _Srcs>
	using _Tuple_implicit_t = enable_if_t<_Tuple_implicit_val<_Dest, _Srcs...>::value, int>;

	
template<bool _Same,
	class _Dest,
	class... _Srcs>
	struct _Tuple_explicit_val0
		: false_type
	{	
	};

template<class... _Dests,
	class... _Srcs>
	struct _Tuple_explicit_val0<true, tuple<_Dests...>, _Srcs...>
		: conjunction<
			is_constructible<_Dests, _Srcs>...,
			negation<conjunction<is_convertible<_Srcs, _Dests>...>>
		>::type
	{	
	};

template<class _Dest,
	class... _Srcs>
	struct _Tuple_explicit_val
		: _Tuple_explicit_val0<tuple_size<_Dest>::value == sizeof...(_Srcs), _Dest, _Srcs...>::type
	{	
	};

template<class _Dest,
	class... _Srcs>
	using _Tuple_explicit_t = enable_if_t<_Tuple_explicit_val<_Dest, _Srcs...>::value, int>;

	
template<class _Myt,
	class... _Other>
	struct _Tuple_convert_copy
	{	
	typedef int type;
	};

template<class _This,
	class _Uty>
	struct _Tuple_convert_copy<tuple<_This>, _Uty>
		: enable_if<!is_same<_This, _Uty>::value
			&& !is_constructible<_This, const tuple<_Uty>&>::value
			&& !is_convertible<const tuple<_Uty>&, _This>::value, int>
	{	
	};

template<class _Myt,
	class... _Other>
	using _Tuple_convert_copy_t = typename _Tuple_convert_copy<_Myt, _Other...>::type;

	
template<class _Myt,
	class... _Other>
	struct _Tuple_convert_move
	{	
	typedef int type;
	};

template<class _This,
	class _Uty>
	struct _Tuple_convert_move<tuple<_This>, _Uty>
		: enable_if<!is_same<_This, _Uty>::value
			&& !is_constructible<_This, tuple<_Uty> >::value
			&& !is_convertible<tuple<_Uty>, _This>::value, int>
	{	
	};

template<class _Myt,
	class... _Other>
	using _Tuple_convert_move_t = typename _Tuple_convert_move<_Myt, _Other...>::type;

	
template<class _Myt,
	class _This2,
	class... _Rest2>
	struct _Tuple_perfect_val
		: true_type
	{	
	};

template<class _Myt,
	class _This2>
	struct _Tuple_perfect_val<_Myt, _This2>
		: negation<is_same<_Myt, remove_const_t<remove_reference_t<_This2>>>>::type
	{	
	};

	
struct _Ignore
	{	
	template<class _Ty>
		void operator=(const _Ty&) const
		{	
		}
	};

constexpr _Ignore ignore{};

		
struct _Tuple_alloc_t
	{	
	};

constexpr _Tuple_alloc_t _Tuple_alloc{};

	
template<class _Ty>
	struct _Tuple_val
	{	
	constexpr _Tuple_val()
		: _Val()
		{	
		}

	template<class _Other>
		constexpr _Tuple_val(_Other&& _Arg)
		: _Val(::std:: forward<_Other>(_Arg))
		{	
		}

	template<class _Other>
		_Tuple_val& operator=(_Other&& _Right)
		{	
		_Val = ::std:: forward<_Other>(_Right);
		return (*this);
		}

	template<class _Alloc,
		class... _Other>
		_Tuple_val(const _Alloc&,
			typename enable_if<!uses_allocator<_Ty, _Alloc>::value,
				_Tuple_alloc_t>::type, _Other&&... _Arg)
		: _Val(::std:: forward<_Other>(_Arg)...)
		{	
		}

	template<class _Alloc,
		class... _Other>
		_Tuple_val(const _Alloc& _Al,
			typename enable_if<uses_allocator<_Ty, _Alloc>::value
				&& is_constructible<_Ty,
					allocator_arg_t, _Alloc>::value,
				_Tuple_alloc_t>::type, _Other&&... _Arg)
		: _Val(allocator_arg, _Al, ::std:: forward<_Other>(_Arg)...)
		{	
		}

	template<class _Alloc,
		class... _Other>
		_Tuple_val(const _Alloc& _Al,
			typename enable_if<uses_allocator<_Ty, _Alloc>::value
				&& !is_constructible<_Ty,
					allocator_arg_t, _Alloc>::value,
				_Tuple_alloc_t>::type, _Other&&... _Arg)
		: _Val(::std:: forward<_Other>(_Arg)..., _Al)
		{	
		}

	_Ty _Val;
	};

	
struct _Exact_args_t
	{	
	};

struct _Unpack_tuple_t
	{	
	};

struct _Alloc_exact_args_t
	{	
	};

struct _Alloc_unpack_tuple_t
	{	
	};

template<class... _Types>
	class tuple;

template<>
	class tuple<>
	{	
public:
	typedef tuple<> _Myt;

	constexpr tuple() noexcept
		{	
		}

	template<class _Alloc>
		tuple(allocator_arg_t, const _Alloc&) noexcept
		{	
		}

	constexpr tuple(const tuple&) noexcept
		{	
		}

	template<class _Alloc>
		tuple(allocator_arg_t, const _Alloc&, const _Myt&) noexcept
		{	
		}

	template<class _Tag,
		enable_if_t<is_same<_Tag, _Exact_args_t>::value, int> = 0>
		constexpr tuple(_Tag) noexcept
		{	
		}

	template<class _Tag,
		enable_if_t<is_same<_Tag, _Unpack_tuple_t>::value, int> = 0>
		constexpr tuple(_Tag, const _Myt&) noexcept
		{	
		}

	template<class _Tag,
		class _Alloc,
		enable_if_t<is_same<_Tag, _Alloc_exact_args_t>::value, int> = 0>
		tuple(_Tag, const _Alloc&) noexcept
		{	
		}

	void swap(_Myt&) noexcept
		{	
		}

	constexpr bool _Equals(const _Myt&) const noexcept
		{	
		return (true);
		}

	constexpr bool _Less(const _Myt&) const noexcept
		{	
		return (false);
		}
	};

template<class _This,
	class... _Rest>
	class tuple<_This, _Rest...>
		: private tuple<_Rest...>
	{	
public:
	typedef _This _This_type;
	typedef tuple<_This, _Rest...> _Myt;
	typedef tuple<_Rest...> _Mybase;
	static constexpr size_t _Mysize = 1 + sizeof...(_Rest);

	template<class _Tag,
		class _This2,
		class... _Rest2,
		enable_if_t<is_same<_Tag, _Exact_args_t>::value, int> = 0>
		constexpr tuple(_Tag, _This2&& _This_arg, _Rest2&&... _Rest_arg)
		: _Mybase(_Exact_args_t{}, ::std:: forward<_Rest2>(_Rest_arg)...),
			_Myfirst(::std:: forward<_This2>(_This_arg))
		{	
		}




















	template<class _Tag,
		class... _Other,
		enable_if_t<is_same<_Tag, _Unpack_tuple_t>::value, int> = 0>
		constexpr tuple(_Tag, const tuple<_Other...>& _Right)
		: _Mybase(_Unpack_tuple_t{}, _Right._Get_rest()),
			_Myfirst(_Right._Myfirst._Val)
		{	
		}

	template<class _Tag,
		class... _Other,
		enable_if_t<is_same<_Tag, _Unpack_tuple_t>::value, int> = 0>
		constexpr tuple(_Tag, tuple<_Other...>&& _Right)
		: _Mybase(_Unpack_tuple_t{}, (typename tuple<_Other...>::_Mybase&&) _Right),
			_Myfirst(::std:: forward<typename tuple<_Other...>::_This_type>(_Right._Myfirst._Val))
		{	
		}

	template<class _Tag,
		class _Alloc,
		class _This2,
		class... _Rest2,
		enable_if_t<is_same<_Tag, _Alloc_exact_args_t>::value, int> = 0>
		tuple(_Tag, const _Alloc& _Al, _This2&& _This_arg, _Rest2&&... _Rest_arg)
		: _Mybase(_Alloc_exact_args_t{}, _Al, ::std:: forward<_Rest2>(_Rest_arg)...),
			_Myfirst(_Al, _Tuple_alloc, ::std:: forward<_This2>(_This_arg))
		{	
		}

	template<class _Tag,
		class _Alloc,
		class _Tpl,
		size_t... _Indices,
		enable_if_t<is_same<_Tag, _Alloc_unpack_tuple_t>::value, int> = 0> inline
		tuple(_Tag, const _Alloc& _Al, _Tpl&& _Right, integer_sequence<size_t, _Indices...>);

	template<class _Tag,
		class _Alloc,
		class _Tpl,
		enable_if_t<is_same<_Tag, _Alloc_unpack_tuple_t>::value, int> = 0>
		tuple(_Tag, const _Alloc& _Al, _Tpl&& _Right)
		: tuple(_Alloc_unpack_tuple_t{}, _Al, ::std:: forward<_Tpl>(_Right),
			make_integer_sequence<size_t, tuple_size<remove_reference_t<_Tpl>>::value>{})
		{	
		}

	template<class _This2 = _This,
		class = enable_if_t<conjunction<is_default_constructible<_This2>,
										is_default_constructible<_Rest>...>::value> >
		constexpr tuple()
		: _Mybase(), _Myfirst()
		{	
		}

	template<class... _Other,
		_Tuple_implicit_t<_Myt, const _Other&...> = 0,
		_Tuple_convert_copy_t<_Myt, _Other...> = 0>
		constexpr tuple(const tuple<_Other...>& _Right)


		: _Mybase(_Unpack_tuple_t{}, _Right._Get_rest()),
			_Myfirst(_Right._Myfirst._Val)
		{	
		}

	template<class... _Other,
		_Tuple_explicit_t<_Myt, const _Other&...> = 0,
		_Tuple_convert_copy_t<_Myt, _Other...> = 0>
		constexpr explicit tuple(const tuple<_Other...>& _Right)


		: _Mybase(_Unpack_tuple_t{}, _Right._Get_rest()),
			_Myfirst(_Right._Myfirst._Val)
		{	
		}

	template<class _Alloc,
		class... _Other,
		_Tuple_implicit_t<_Myt, const _Other&...> = 0,
		_Tuple_convert_copy_t<_Myt, _Other...> = 0>
		tuple(allocator_arg_t, const _Alloc& _Al,
			const tuple<_Other...>& _Right)
		: tuple(_Alloc_unpack_tuple_t{}, _Al, _Right)
		{	
		}

	template<class _Alloc,
		class... _Other,
		_Tuple_explicit_t<_Myt, const _Other&...> = 0,
		_Tuple_convert_copy_t<_Myt, _Other...> = 0>
		explicit tuple(allocator_arg_t, const _Alloc& _Al,
			const tuple<_Other...>& _Right)
		: tuple(_Alloc_unpack_tuple_t{}, _Al, _Right)
		{	
		}

	template<class _This2 = _This,
		_Tuple_implicit_t<_Myt, const _This2&, const _Rest&...> = 0>
		constexpr tuple(const _This& _This_arg, const _Rest&... _Rest_arg)


		: _Mybase(_Exact_args_t{}, _Rest_arg...), _Myfirst(_This_arg)
		{	
		}

	template<class _This2 = _This,
		_Tuple_explicit_t<_Myt, const _This2&, const _Rest&...> = 0>
		constexpr explicit tuple(const _This& _This_arg, const _Rest&... _Rest_arg)


		: _Mybase(_Exact_args_t{}, _Rest_arg...), _Myfirst(_This_arg)
		{	
		}

	template<class _Alloc,
		class _This2 = _This,
		_Tuple_implicit_t<_Myt, const _This2&, const _Rest&...> = 0>
		tuple(allocator_arg_t, const _Alloc& _Al,
			const _This& _This_arg, const _Rest&... _Rest_arg)
		: tuple(_Alloc_exact_args_t{}, _Al, _This_arg, _Rest_arg...)
		{	
		}

	template<class _Alloc,
		class _This2 = _This,
		_Tuple_explicit_t<_Myt, const _This2&, const _Rest&...> = 0>
		explicit tuple(allocator_arg_t, const _Alloc& _Al,
			const _This& _This_arg, const _Rest&... _Rest_arg)
		: tuple(_Alloc_exact_args_t{}, _Al, _This_arg, _Rest_arg...)
		{	
		}

	template<class _This2,
		class... _Rest2,
		enable_if_t<conjunction<
			_Tuple_perfect_val<_Myt, _This2, _Rest2...>,
			_Tuple_implicit_val<_Myt, _This2, _Rest2...>
		>::value, int> = 0>
		constexpr tuple(_This2&& _This_arg, _Rest2&&... _Rest_arg)


		: _Mybase(_Exact_args_t{}, ::std:: forward<_Rest2>(_Rest_arg)...),
			_Myfirst(::std:: forward<_This2>(_This_arg))
		{	
		}

	template<class _This2,
		class... _Rest2,
		enable_if_t<conjunction<
			_Tuple_perfect_val<_Myt, _This2, _Rest2...>,
			_Tuple_explicit_val<_Myt, _This2, _Rest2...>
		>::value, int> = 0>
		constexpr explicit tuple(_This2&& _This_arg, _Rest2&&... _Rest_arg)


		: _Mybase(_Exact_args_t{}, ::std:: forward<_Rest2>(_Rest_arg)...),
			_Myfirst(::std:: forward<_This2>(_This_arg))
		{	
		}

	template<class _Alloc,
		class _This2,
		class... _Rest2,
		enable_if_t<conjunction<
			_Tuple_perfect_val<_Myt, _This2, _Rest2...>,
			_Tuple_implicit_val<_Myt, _This2, _Rest2...>
		>::value, int> = 0>
		tuple(allocator_arg_t, const _Alloc& _Al,
			_This2&& _This_arg, _Rest2&&... _Rest_arg)
		: tuple(_Alloc_exact_args_t{}, _Al, ::std:: forward<_This2>(_This_arg), ::std:: forward<_Rest2>(_Rest_arg)...)
		{	
		}

	template<class _Alloc,
		class _This2,
		class... _Rest2,
		enable_if_t<conjunction<
			_Tuple_perfect_val<_Myt, _This2, _Rest2...>,
			_Tuple_explicit_val<_Myt, _This2, _Rest2...>
		>::value, int> = 0>
		explicit tuple(allocator_arg_t, const _Alloc& _Al,
			_This2&& _This_arg, _Rest2&&... _Rest_arg)
		: tuple(_Alloc_exact_args_t{}, _Al, ::std:: forward<_This2>(_This_arg), ::std:: forward<_Rest2>(_Rest_arg)...)
		{	
		}

	template<class... _Other,
		_Tuple_implicit_t<_Myt, _Other...> = 0,
		_Tuple_convert_move_t<_Myt, _Other...> = 0>
		constexpr tuple(tuple<_Other...>&& _Right)


		: _Mybase(_Unpack_tuple_t{}, (typename tuple<_Other...>::_Mybase&&) _Right),
			_Myfirst(::std:: forward<typename tuple<_Other...>::_This_type>(_Right._Myfirst._Val))
		{	
		}

	template<class... _Other,
		_Tuple_explicit_t<_Myt, _Other...> = 0,
		_Tuple_convert_move_t<_Myt, _Other...> = 0>
		constexpr explicit tuple(tuple<_Other...>&& _Right)


		: _Mybase(_Unpack_tuple_t{}, (typename tuple<_Other...>::_Mybase&&) _Right),
			_Myfirst(::std:: forward<typename tuple<_Other...>::_This_type>(_Right._Myfirst._Val))
		{	
		}

	template<class _Alloc,
		class... _Other,
		_Tuple_implicit_t<_Myt, _Other...> = 0,
		_Tuple_convert_move_t<_Myt, _Other...> = 0>
		tuple(allocator_arg_t, const _Alloc& _Al,
			tuple<_Other...>&& _Right)
		: tuple(_Alloc_unpack_tuple_t{}, _Al, ::std:: move(_Right))
		{	
		}

	template<class _Alloc,
		class... _Other,
		_Tuple_explicit_t<_Myt, _Other...> = 0,
		_Tuple_convert_move_t<_Myt, _Other...> = 0>
		explicit tuple(allocator_arg_t, const _Alloc& _Al,
			tuple<_Other...>&& _Right)
		: tuple(_Alloc_unpack_tuple_t{}, _Al, ::std:: move(_Right))
		{	
		}

	template<class... _Other>
		_Myt& operator=(const tuple<_Other...>& _Right)
		{	
		_Myfirst._Val = _Right._Myfirst._Val;
		_Get_rest() = _Right._Get_rest();
		return (*this);
		}

	template<class... _Other>
		_Myt& operator=(tuple<_Other...>&& _Right)
		{	
		_Myfirst._Val = ::std:: forward<typename tuple<_Other...>::_This_type>
			(_Right._Myfirst._Val);
		_Get_rest() = ::std:: forward<typename tuple<_Other...>::_Mybase>
			(_Right._Get_rest());
		return (*this);
		}

	template<class... _Other>
		constexpr bool _Equals(const tuple<_Other...>& _Right) const
		{	
		static_assert(_Mysize == sizeof...(_Other),
			"comparing tuple to object with different size");
		return (_Myfirst._Val == _Right._Myfirst._Val
			&& _Mybase::_Equals(_Right._Get_rest()));
		}

	template<class... _Other>
		constexpr bool _Less(const tuple<_Other...>& _Right) const
		{	
		static_assert(_Mysize == sizeof...(_Other),
			"comparing tuple to object with different size");
		return (_Myfirst._Val < _Right._Myfirst._Val
			|| (!(_Right._Myfirst._Val < _Myfirst._Val)
				&& _Mybase::_Less(_Right._Get_rest())));
		}

	template<class _Alloc,
		class _This2 = _This,
		class = enable_if_t<conjunction<is_default_constructible<_This2>,
										is_default_constructible<_Rest>...>::value> >
		tuple(allocator_arg_t, const _Alloc& _Al)
		: _Mybase(allocator_arg, _Al), _Myfirst(_Al, _Tuple_alloc)
		{	
		}

	template<class _Alloc>
		tuple(allocator_arg_t, const _Alloc& _Al,
			const _Myt& _Right)
		: tuple(_Alloc_unpack_tuple_t{}, _Al, _Right)
		{	
		}

	tuple(const _Myt&) = default;
	tuple(_Myt&&) = default;

	template<class _First,
		class _Second,
		_Tuple_implicit_t<_Myt, const _First&, const _Second&> = 0>
		constexpr tuple(const pair<_First, _Second>& _Right)


		: _Mybase(_Exact_args_t{}, _Right.second), _Myfirst(_Right.first)
		{	
		}

	template<class _First,
		class _Second,
		_Tuple_explicit_t<_Myt, const _First&, const _Second&> = 0>
		constexpr explicit tuple(const pair<_First, _Second>& _Right)


		: _Mybase(_Exact_args_t{}, _Right.second), _Myfirst(_Right.first)
		{	
		}

	template<class _Alloc,
		class _First,
		class _Second,
		_Tuple_implicit_t<_Myt, const _First&, const _Second&> = 0>
		tuple(allocator_arg_t, const _Alloc& _Al,
			const pair<_First, _Second>& _Right)
		: tuple(_Alloc_unpack_tuple_t{}, _Al, _Right)
		{	
		}

	template<class _Alloc,
		class _First,
		class _Second,
		_Tuple_explicit_t<_Myt, const _First&, const _Second&> = 0>
		explicit tuple(allocator_arg_t, const _Alloc& _Al,
			const pair<_First, _Second>& _Right)
		: tuple(_Alloc_unpack_tuple_t{}, _Al, _Right)
		{	
		}

	_Myt& operator=(const _Myt& _Right)
		{	
		_Myfirst._Val = _Right._Myfirst._Val;
		_Get_rest() = _Right._Get_rest();
		return (*this);
		}

	template<class _First,
		class _Second>
		_Myt& operator=(const pair<_First, _Second>& _Right)
		{	
		static_assert(_Mysize == 2,
			"assigning to tuple from object with different size");
		_Myfirst._Val = _Right.first;
		_Get_rest()._Myfirst._Val = _Right.second;
		return (*this);
		}

	template<class _Alloc>
		tuple(allocator_arg_t, const _Alloc& _Al,
			_Myt&& _Right)
		: tuple(_Alloc_unpack_tuple_t{}, _Al, ::std:: move(_Right))
		{	
		}

	template<class _First,
		class _Second,
		_Tuple_implicit_t<_Myt, _First, _Second> = 0>
		constexpr tuple(pair<_First, _Second>&& _Right)


		: _Mybase(_Exact_args_t{}, ::std:: forward<_Second>(_Right.second)),
			_Myfirst(::std:: forward<_First>(_Right.first))
		{	
		}

	template<class _First,
		class _Second,
		_Tuple_explicit_t<_Myt, _First, _Second> = 0>
		constexpr explicit tuple(pair<_First, _Second>&& _Right)


		: _Mybase(_Exact_args_t{}, ::std:: forward<_Second>(_Right.second)),
			_Myfirst(::std:: forward<_First>(_Right.first))
		{	
		}

	template<class _Alloc,
		class _First,
		class _Second,
		_Tuple_implicit_t<_Myt, _First, _Second> = 0>
		tuple(allocator_arg_t, const _Alloc& _Al,
			pair<_First, _Second>&& _Right)
		: tuple(_Alloc_unpack_tuple_t{}, _Al, ::std:: move(_Right))
		{	
		}

	template<class _Alloc,
		class _First,
		class _Second,
		_Tuple_explicit_t<_Myt, _First, _Second> = 0>
		explicit tuple(allocator_arg_t, const _Alloc& _Al,
			pair<_First, _Second>&& _Right)
		: tuple(_Alloc_unpack_tuple_t{}, _Al, ::std:: move(_Right))
		{	
		}

	_Myt& operator=(_Myt&& _Right)
		noexcept(is_nothrow_move_assignable<_This>::value && is_nothrow_move_assignable<_Mybase>::value)
		{	
		_Myfirst._Val = ::std:: forward<_This>(_Right._Myfirst._Val);
		_Get_rest() = ::std:: forward<_Mybase>(_Right._Get_rest());
		return (*this);
		}

	template<class _First,
		class _Second>
		_Myt& operator=(pair<_First, _Second>&& _Right)
		{	
		static_assert(_Mysize == 2,
			"assigning to tuple from object with different size");
		_Myfirst._Val = ::std:: forward<_First>(_Right.first);
		_Get_rest()._Myfirst._Val = ::std:: forward<_Second>(_Right.second);
		return (*this);
		}

	_Mybase& _Get_rest() noexcept
		{	
		return (*this);
		}

	constexpr const _Mybase& _Get_rest() const noexcept
		{	
		return (*this);
		}

	_Tuple_val<_This> _Myfirst;	

	void swap(tuple& _Right)
		noexcept((conjunction<_Is_nothrow_swappable<_This>, _Is_nothrow_swappable<_Rest>...>::value))
		{	
		_Swap_adl(_Myfirst._Val, _Right._Myfirst._Val);
		_Mybase::swap(_Right._Get_rest());
		}
	};


	

template<class... _Types1,
	class... _Types2> inline
	constexpr bool operator==(const tuple<_Types1...>& _Left,
		const tuple<_Types2...>& _Right)
	{	
	return (_Left._Equals(_Right));
	}

template<class... _Types1,
	class... _Types2> inline
	constexpr bool operator!=(const tuple<_Types1...>& _Left,
		const tuple<_Types2...>& _Right)
	{	
	return (!(_Left == _Right));
	}

template<class... _Types1,
	class... _Types2> inline
	constexpr bool operator<(const tuple<_Types1...>& _Left,
		const tuple<_Types2...>& _Right)
	{	
	return (_Left._Less(_Right));
	}

template<class... _Types1,
	class... _Types2> inline
	constexpr bool operator>=(const tuple<_Types1...>& _Left,
		const tuple<_Types2...>& _Right)
	{	
	return (!(_Left < _Right));
	}

template<class... _Types1,
	class... _Types2> inline
	constexpr bool operator>(const tuple<_Types1...>& _Left,
		const tuple<_Types2...>& _Right)
	{	
	return (_Right < _Left);
	}

template<class... _Types1,
	class... _Types2> inline
	constexpr bool operator<=(const tuple<_Types1...>& _Left,
		const tuple<_Types2...>& _Right)
	{	
	return (!(_Right < _Left));
	}

template<class... _Types,
	class = enable_if_t<conjunction<_Is_swappable<_Types>...>::value>> inline
	void swap(tuple<_Types...>& _Left,
		tuple<_Types...>& _Right)
			noexcept(noexcept(_Left.swap(_Right)))
	{	
	return (_Left.swap(_Right));
	}


	
template<class _Ty,
	class _Tuple>
	struct _Tuple_element;

template<class _This,
	class... _Rest>
	struct _Tuple_element<_This, tuple<_This, _Rest...> >
	{	
	typedef int _Check_type;
	static_assert(is_void<typename _Tuple_element<_This,
		tuple<_Rest...> >::_Check_type>::value,
		"duplicate type T in get<T>(tuple)");

	typedef _This type;
	typedef tuple<_This, _Rest...> _Ttype;
	};

template<class _Ty,
	class _This,
	class... _Rest>
	struct _Tuple_element<_Ty, tuple<_This, _Rest...> >
		: public _Tuple_element<_Ty, tuple<_Rest...> >
	{	
	};

template<class _Ty>
	struct _Tuple_element<_Ty, tuple<> >
	{	
	typedef void _Check_type;	
	};

template<class _Ty,
	class _Tuple>
	struct _Tuple_element<_Ty, const _Tuple>
		: public _Tuple_element<_Ty, _Tuple>
	{	
	typedef _Tuple_element<_Ty, _Tuple> _Mybase;
	typedef typename add_const<typename _Mybase::type>::type type;
	};

template<class _Ty,
	class _Tuple>
	struct _Tuple_element<_Ty, volatile _Tuple>
		: public _Tuple_element<_Ty, _Tuple>
	{	
	typedef _Tuple_element<_Ty, _Tuple> _Mybase;
	typedef typename add_volatile<typename _Mybase::type>::type type;
	};

template<class _Ty,
	class _Tuple>
	struct _Tuple_element<_Ty, const volatile _Tuple>
		: public _Tuple_element<_Ty, _Tuple>
	{	
	typedef _Tuple_element<_Ty, _Tuple> _Mybase;
	typedef typename add_cv<typename _Mybase::type>::type type;
	};

	
template<size_t _Index,
	class... _Types> inline
	constexpr typename tuple_element<_Index, tuple<_Types...> >::type&
		get(tuple<_Types...>& _Tuple) noexcept
	{	
	typedef typename tuple_element<_Index, tuple<_Types...> >::_Ttype
		_Ttype;
	return (((_Ttype&)_Tuple)._Myfirst._Val);
	}

template<size_t _Index,
	class... _Types> inline
	constexpr const typename tuple_element<_Index, tuple<_Types...> >::type&
		get(const tuple<_Types...>& _Tuple) noexcept
	{	
	typedef typename tuple_element<_Index, tuple<_Types...> >::_Ttype
		_Ttype;
	return (((_Ttype&)_Tuple)._Myfirst._Val);
	}

template<size_t _Index,
	class... _Types> inline
	constexpr typename tuple_element<_Index, tuple<_Types...> >::type&&
		get(tuple<_Types...>&& _Tuple) noexcept
	{	
	typedef typename tuple_element<_Index, tuple<_Types...> >::_Ttype
		_Ttype;
	typedef typename tuple_element<_Index, tuple<_Types...> >::type&&
		_RRtype;
	return (::std:: forward<_RRtype>(((_Ttype&)_Tuple)._Myfirst._Val));
	}

	
template<class _Ty,
	class... _Types> inline
	constexpr _Ty& get(tuple<_Types...>& _Tuple) noexcept
	{	
	typedef typename _Tuple_element<_Ty, tuple<_Types...> >::_Ttype _Ttype;
	return (((_Ttype&)_Tuple)._Myfirst._Val);
	}

template<class _Ty,
	class... _Types> inline
	constexpr const _Ty& get(const tuple<_Types...>& _Tuple) noexcept
	{	
	typedef typename _Tuple_element<_Ty, tuple<_Types...> >::_Ttype _Ttype;
	return (((_Ttype&)_Tuple)._Myfirst._Val);
	}

template<class _Ty,
	class... _Types> inline
	constexpr _Ty&& get(tuple<_Types...>&& _Tuple) noexcept
	{	
	typedef typename _Tuple_element<_Ty, tuple<_Types...> >::_Ttype _Ttype;
	return (::std:: forward<_Ty&&>(((_Ttype&)_Tuple)._Myfirst._Val));
	}

	













template<class _This,
	class... _Rest>
	template<class _Tag,
		class _Alloc,
		class _Tpl,
		size_t... _Indices,
		enable_if_t<is_same<_Tag, _Alloc_unpack_tuple_t>::value, int>> inline
		tuple<_This, _Rest...>::tuple(_Tag, const _Alloc& _Al, _Tpl&& _Right, integer_sequence<size_t, _Indices...>)
		: tuple(_Alloc_exact_args_t{}, _Al, ::std:: get<_Indices>(::std:: forward<_Tpl>(_Right))...)
		{	
		}

	
template<class... _Types> inline
	constexpr tuple<typename _Unrefwrap<_Types>::type...>
		make_tuple(_Types&&... _Args)
	{	
	typedef tuple<typename _Unrefwrap<_Types>::type...> _Ttype;
	return (_Ttype(::std:: forward<_Types>(_Args)...));
	}

	
template<class... _Types> inline
	constexpr tuple<_Types&...>
		tie(_Types&... _Args) noexcept
	{	
	typedef tuple<_Types&...> _Ttype;
	return (_Ttype(_Args...));
	}


	

template<class... _Types> inline
	constexpr tuple<_Types&&...>
		forward_as_tuple(_Types&&... _Args) noexcept
	{	
	return (tuple<_Types&&...>(::std:: forward<_Types>(_Args)...));
	}


	
template<class _Seq_type1,
	class _Seq_type2>
	struct _Cat_sequences;

template<size_t... _Indexes1,
	size_t... _Indexes2>
	struct _Cat_sequences<integer_sequence<size_t, _Indexes1...>,
		integer_sequence<size_t, _Indexes2...> >
	{	
	typedef integer_sequence<size_t, _Indexes1..., _Indexes2...> type;
	};

	
template<class _Ty,
	size_t _Size>
	class array;

template<size_t _Idx,
	class _Ty,
	size_t _Size>
	constexpr _Ty& get(array<_Ty, _Size>& _Arr) noexcept;

template<size_t _Idx,
	class _Ty,
	size_t _Size>
	constexpr const _Ty& get(const array<_Ty, _Size>& _Arr) noexcept;

template<size_t _Idx,
	class _Ty,
	size_t _Size>
	constexpr _Ty&& get(array<_Ty, _Size>&& _Arr) noexcept;

	
template<class _Ty,
	class... _For_array>
	struct _View_as_tuple
	{	
	static_assert(_Always_false<_Ty>::value,
		"Unsupported tuple_cat arguments.");
	};

template<class... _Types>
	struct _View_as_tuple<tuple<_Types...> >
	{	
	typedef tuple<_Types...> type;
	};

template<class _Ty1,
	class _Ty2>
	struct _View_as_tuple<pair<_Ty1, _Ty2> >
	{	
	typedef tuple<_Ty1, _Ty2> type;
	};

template<class _Ty,
	class... _Types>
	struct _View_as_tuple<array<_Ty, 0>, _Types...>
	{	
	typedef tuple<_Types...> type;
	};

template<class _Ty,
	size_t _Size,
	class... _Types>
	struct _View_as_tuple<array<_Ty, _Size>, _Types...>
		: _View_as_tuple<array<_Ty, _Size - 1>, _Ty, _Types...>
	{	
	};

	
template<size_t _Nx,
	class _Ty>
	struct _Repeat_for
		: integral_constant<size_t, _Nx>
	{	
	};

	
template<class _Ret,
	class _Kx_arg,
	class _Ix_arg,
	size_t _Ix_next,
	class... _Tuples>
	struct _Tuple_cat2
	{	
	static_assert(sizeof...(_Tuples) == 0,
		"Unsupported tuple_cat arguments.");
	typedef _Ret type;
	typedef _Kx_arg _Kx_arg_seq;
	typedef _Ix_arg _Ix_arg_seq;
	};

template<class... _Types1,
	class _Kx_arg,
	size_t... _Ix,
	size_t _Ix_next,
	class... _Types2,
	class... _Rest>
	struct _Tuple_cat2<tuple<_Types1...>, _Kx_arg,
		integer_sequence<size_t, _Ix...>, _Ix_next,
		tuple<_Types2...>, _Rest...>
		: _Tuple_cat2<
			tuple<_Types1..., _Types2...>,
			typename _Cat_sequences<_Kx_arg,
				make_integer_sequence<size_t, sizeof...(_Types2)> >::type,
			integer_sequence<size_t, _Ix...,
				_Repeat_for<_Ix_next, _Types2>::value...>,
			_Ix_next + 1,
			_Rest...>
	{	
	};

template<class... _Tuples>
	struct _Tuple_cat1
		: _Tuple_cat2<tuple<>, integer_sequence<size_t>,
				integer_sequence<size_t>, 0,
			typename _View_as_tuple<typename decay<_Tuples>::type>::type...>
	{	
	};

template<class _Ret,
	size_t... _Kx,
	size_t... _Ix,
	class _Ty> inline
	constexpr _Ret _Tuple_cat(integer_sequence<size_t, _Kx...>,
		integer_sequence<size_t, _Ix...>, _Ty&& _Arg)
	{	
	return (_Ret(::std:: get<_Kx>(::std:: get<_Ix>(::std:: forward<_Ty>(_Arg)))...));
	}

template<class _Ret,
	class _Ty> inline
	constexpr _Ret _Tuple_cat(integer_sequence<size_t>,
		integer_sequence<size_t>, _Ty&&)
	{	
	return (_Ret());
	}

template<class... _Tuples> inline
	constexpr typename _Tuple_cat1<_Tuples...>::type
		tuple_cat(_Tuples&&... _Tpls)
	{	
	typedef _Tuple_cat1<_Tuples...> _Cat1;
	return (_Tuple_cat<typename _Cat1::type>(
		typename _Cat1::_Kx_arg_seq(), typename _Cat1::_Ix_arg_seq(),
		::std:: forward_as_tuple(::std:: forward<_Tuples>(_Tpls)...)));
	}


	
template<class _Tpl,
	class _Fx,
	size_t... _Indices> inline
	void _For_each_tuple_element_impl(_Tpl&& _Tuple,
		_Fx _Func, integer_sequence<size_t, _Indices...>)
	{	
	int _Ignored[] = { (static_cast<void>(_Func(
		::std:: get<_Indices>(::std:: forward<_Tpl>(_Tuple))
		)), 0)... };
	(void)_Ignored;
	}

template<class _Tpl,
	class _Fx> inline
	void _For_each_tuple_element(_Tpl&& _Tuple, _Fx _Func)
	{	
	_For_each_tuple_element_impl(
		::std:: forward<_Tpl>(_Tuple),
		_Func,
		make_integer_sequence<size_t,
			tuple_size<remove_reference_t<_Tpl>>::value>()
		);
	}


	
template<class _Ty1,
	class _Ty2>
	template<class _Tuple1,
		class _Tuple2,
		size_t... _Indexes1,
		size_t... _Indexes2> inline
		pair<_Ty1, _Ty2>::pair(_Tuple1& _Val1,
			_Tuple2& _Val2,
			integer_sequence<size_t, _Indexes1...>,
			integer_sequence<size_t, _Indexes2...>)
		: first(::std:: get<_Indexes1>(::std:: move(_Val1))...),
			second(::std:: get<_Indexes2>(::std:: move(_Val2))...)
		{	
		(void) _Val1;	
		(void) _Val2;
		}

	
template<class _Ty1,
	class _Ty2>
	template<class... _Types1,
		class... _Types2> inline
		pair<_Ty1, _Ty2>::pair(piecewise_construct_t,
			tuple<_Types1...> _Val1,
			tuple<_Types2...> _Val2)
		: pair(_Val1, _Val2,
			make_integer_sequence<size_t, sizeof...(_Types1)>(),
			make_integer_sequence<size_t, sizeof...(_Types2)>())
		{	
		}

}

namespace std {
	
template<class... _Types,
	class _Alloc>
	struct uses_allocator<tuple<_Types...>, _Alloc>
		: true_type
	{	
	};

}	


namespace std {
namespace tr1 {	
using ::std:: get;
using ::std:: ignore;
using ::std:: make_tuple;
using ::std:: ref;
using ::std:: tie;
using ::std:: tuple;
}	
}


 
 #pragma warning(pop)
 #pragma pack(pop)












#pragma once






 #pragma pack(push,8)
 #pragma warning(push,3)
 
 

 #pragma warning(disable: 4127)
namespace std {
		
template<class _Mytree,
	class _Base = _Iterator_base0>
	class _Tree_unchecked_const_iterator
		: public _Iterator012<bidirectional_iterator_tag,
			typename _Mytree::value_type,
			typename _Mytree::difference_type,
			typename _Mytree::const_pointer,
			typename _Mytree::const_reference,
			_Base>
	{	
public:
	typedef _Tree_unchecked_const_iterator<_Mytree, _Base> _Myiter;
	typedef bidirectional_iterator_tag iterator_category;

	typedef typename _Mytree::_Nodeptr _Nodeptr;
	typedef typename _Mytree::value_type value_type;
	typedef typename _Mytree::difference_type difference_type;
	typedef typename _Mytree::const_pointer pointer;
	typedef typename _Mytree::const_reference reference;

	_Tree_unchecked_const_iterator()
		: _Ptr()
		{	
		}

	_Tree_unchecked_const_iterator(_Nodeptr _Pnode, const _Mytree *_Plist)
		: _Ptr(_Pnode)
		{	
		this->_Adopt(_Plist);
		}

	reference operator*() const
		{	
		return (_Mytree::_Myval(_Ptr));
		}

	pointer operator->() const
		{	
		return (pointer_traits<pointer>::pointer_to(**this));
		}

	_Myiter& operator++()
		{	
		if (_Mytree::_Isnil(_Ptr))
			;	
		else if (!_Mytree::_Isnil(_Mytree::_Right(_Ptr)))
			_Ptr = _Mytree::_Min(
				_Mytree::_Right(_Ptr));	
		else
			{	
			_Nodeptr _Pnode;
			while (!_Mytree::_Isnil(_Pnode = _Mytree::_Parent(_Ptr))
				&& _Ptr == _Mytree::_Right(_Pnode))
				_Ptr = _Pnode;	
			_Ptr = _Pnode;	
			}
		return (*this);
		}

	_Myiter operator++(int)
		{	
		_Myiter _Tmp = *this;
		++*this;
		return (_Tmp);
		}

	_Myiter& operator--()
		{	
		if (_Mytree::_Isnil(_Ptr))
			_Ptr = _Mytree::_Right(_Ptr);	
		else if (!_Mytree::_Isnil(_Mytree::_Left(_Ptr)))
			_Ptr = _Mytree::_Max(
				_Mytree::_Left(_Ptr));	
		else
			{	
			_Nodeptr _Pnode;
			while (!_Mytree::_Isnil(_Pnode = _Mytree::_Parent(_Ptr))
				&& _Ptr == _Mytree::_Left(_Pnode))
				_Ptr = _Pnode;	
			if (_Mytree::_Isnil(_Ptr))
				;	
			else
				_Ptr = _Pnode;	
			}
		return (*this);
		}

	_Myiter operator--(int)
		{	
		_Myiter _Tmp = *this;
		--*this;
		return (_Tmp);
		}

	bool operator==(const _Myiter& _Right) const
		{	
		return (_Ptr == _Right._Ptr);
		}

	bool operator!=(const _Myiter& _Right) const
		{	
		return (!(*this == _Right));
		}

	_Nodeptr _Mynode() const
		{	
		return (_Ptr);
		}

	_Nodeptr _Ptr;	
	};

	
template<class _Mytree>
	class _Tree_unchecked_iterator
		: public _Tree_unchecked_const_iterator<_Mytree>
	{	
public:
	typedef _Tree_unchecked_iterator<_Mytree> _Myiter;
	typedef _Tree_unchecked_const_iterator<_Mytree> _Mybase;
	typedef bidirectional_iterator_tag iterator_category;

	typedef typename _Mytree::_Nodeptr _Nodeptr;
	typedef typename _Mytree::value_type value_type;
	typedef typename _Mytree::difference_type difference_type;
	typedef typename _Mytree::pointer pointer;
	typedef typename _Mytree::reference reference;

	_Tree_unchecked_iterator()
		{	
		}

	_Tree_unchecked_iterator(_Nodeptr _Pnode, const _Mytree *_Plist)
		: _Mybase(_Pnode, _Plist)
		{	
		}

	reference operator*() const
		{	
		return ((reference)**(_Mybase *)this);
		}

	pointer operator->() const
		{	
		return (pointer_traits<pointer>::pointer_to(**this));
		}

	_Myiter& operator++()
		{	
		++static_cast<_Mybase&>(*this);
		return (*this);
		}

	_Myiter operator++(int)
		{	
		_Myiter _Tmp = *this;
		++*this;
		return (_Tmp);
		}

	_Myiter& operator--()
		{	
		--static_cast<_Mybase&>(*this);
		return (*this);
		}

	_Myiter operator--(int)
		{	
		_Myiter _Tmp = *this;
		--*this;
		return (_Tmp);
		}
	};

	
template<class _Mytree>
	class _Tree_const_iterator
		: public _Tree_unchecked_const_iterator<_Mytree, _Iterator_base>
	{	
public:
	typedef _Tree_const_iterator<_Mytree> _Myiter;
	typedef _Tree_unchecked_const_iterator<_Mytree, _Iterator_base> _Mybase;
	typedef bidirectional_iterator_tag iterator_category;

	typedef typename _Mytree::_Nodeptr _Nodeptr;
	typedef typename _Mytree::value_type value_type;
	typedef typename _Mytree::difference_type difference_type;
	typedef typename _Mytree::const_pointer pointer;
	typedef typename _Mytree::const_reference reference;

	_Tree_const_iterator()
		: _Mybase()
		{	
		}

	_Tree_const_iterator(_Nodeptr _Pnode, const _Mytree *_Plist)
		: _Mybase(_Pnode, _Plist)
		{	
		}

	typedef _Tree_unchecked_const_iterator<_Mytree> _Unchecked_type;

	_Myiter& _Rechecked(_Unchecked_type _Right)
		{	
		this->_Ptr = _Right._Ptr;
		return (*this);
		}

	_Unchecked_type _Unchecked() const
		{	
		return (_Unchecked_type(this->_Ptr, (_Mytree *)this->_Getcont()));
		}

	reference operator*() const
		{	
 
















		return (_Mytree::_Myval(this->_Ptr));
		}

	pointer operator->() const
		{	
		return (pointer_traits<pointer>::pointer_to(**this));
		}

	_Myiter& operator++()
		{	
 













		++static_cast<_Mybase&>(*this);
		return (*this);
		}

	_Myiter operator++(int)
		{	
		_Myiter _Tmp = *this;
		++*this;
		return (_Tmp);
		}

	_Myiter& operator--()
		{	
 























		--static_cast<_Mybase&>(*this);
 

		return (*this);
		}

	_Myiter operator--(int)
		{	
		_Myiter _Tmp = *this;
		--*this;
		return (_Tmp);
		}

	bool operator==(const _Myiter& _Right) const
		{	
 










		return (this->_Ptr == _Right._Ptr);
		}

	bool operator!=(const _Myiter& _Right) const
		{	
		return (!(*this == _Right));
		}
	};

template<class _Mytree> inline
	typename _Tree_const_iterator<_Mytree>::_Unchecked_type
		_Unchecked(_Tree_const_iterator<_Mytree> _Iter)
	{	
	return (_Iter._Unchecked());
	}

template<class _Mytree> inline
	_Tree_const_iterator<_Mytree>&
		_Rechecked(_Tree_const_iterator<_Mytree>& _Iter,
			typename _Tree_const_iterator<_Mytree>
				::_Unchecked_type _Right)
	{	
	return (_Iter._Rechecked(_Right));
	}

	
template<class _Mytree>
	class _Tree_iterator
		: public _Tree_const_iterator<_Mytree>
	{	
public:
	typedef _Tree_iterator<_Mytree> _Myiter;
	typedef _Tree_const_iterator<_Mytree> _Mybase;
	typedef bidirectional_iterator_tag iterator_category;

	typedef typename _Mytree::_Nodeptr _Nodeptr;
	typedef typename _Mytree::value_type value_type;
	typedef typename _Mytree::difference_type difference_type;

	typedef typename _Mytree::pointer pointer;
	typedef typename _Mytree::reference reference;

	_Tree_iterator()
		{	
		}

	_Tree_iterator(_Nodeptr _Pnode, const _Mytree *_Plist)
		: _Mybase(_Pnode, _Plist)
		{	
		}

	typedef _Tree_unchecked_iterator<_Mytree> _Unchecked_type;

	_Myiter& _Rechecked(_Unchecked_type _Right)
		{	
		this->_Ptr = _Right._Ptr;
		return (*this);
		}

	_Unchecked_type _Unchecked() const
		{	
		return (_Unchecked_type(this->_Ptr, static_cast<const _Mytree *>(this->_Getcont())));
		}

	reference operator*() const
		{	
		return ((reference)**(_Mybase *)this);
		}

	pointer operator->() const
		{	
		return (pointer_traits<pointer>::pointer_to(**this));
		}

	_Myiter& operator++()
		{	
		++static_cast<_Mybase&>(*this);
		return (*this);
		}

	_Myiter operator++(int)
		{	
		_Myiter _Tmp = *this;
		++*this;
		return (_Tmp);
		}

	_Myiter& operator--()
		{	
		--static_cast<_Mybase&>(*this);
		return (*this);
		}

	_Myiter operator--(int)
		{	
		_Myiter _Tmp = *this;
		--*this;
		return (_Tmp);
		}
	};

template<class _Mytree> inline
	typename _Tree_iterator<_Mytree>::_Unchecked_type
		_Unchecked(_Tree_iterator<_Mytree> _Iter)
	{	
	return (_Iter._Unchecked());
	}

template<class _Mytree> inline
	_Tree_iterator<_Mytree>&
		_Rechecked(_Tree_iterator<_Mytree>& _Iter,
			typename _Tree_iterator<_Mytree>
				::_Unchecked_type _Right)
	{	
	return (_Iter._Rechecked(_Right));
	}

		
template<class _Value_type,
	class _Size_type,
	class _Difference_type,
	class _Pointer,
	class _Const_pointer,
	class _Reference,
	class _Const_reference,
	class _Nodeptr_type>
	struct _Tree_iter_types
	{	
	typedef _Value_type value_type;
	typedef _Size_type size_type;
	typedef _Difference_type difference_type;
	typedef _Pointer pointer;
	typedef _Const_pointer const_pointer;
	typedef _Reference reference;
	typedef _Const_reference const_reference;
	typedef _Nodeptr_type _Nodeptr;
	};

template<class _Value_type,
	class _Voidptr>
	struct _Tree_node
		{	
		_Voidptr _Left;	
		_Voidptr _Parent;	
		_Voidptr _Right;	
		char _Color;	
		char _Isnil;	
		_Value_type _Myval;	

	private:
		_Tree_node& operator=(const _Tree_node&);
		};

template<class _Value_type>
	struct _Tree_node<_Value_type, void *>
		{	
		typedef _Tree_node<_Value_type, void *> *_Nodeptr;
		_Nodeptr _Left;	
		_Nodeptr _Parent;	
		_Nodeptr _Right;	
		char _Color;	
		char _Isnil;	
		_Value_type _Myval;	

	private:
		_Tree_node& operator=(const _Tree_node&);
		};

template<class _Ty>
	struct _Tree_simple_types
		: public _Simple_types<_Ty>
	{	
	typedef _Tree_node<_Ty, void *> _Node;
	typedef _Node *_Nodeptr;
	};

template<class _Ty,
	class _Alloc0>
	struct _Tree_base_types
	{	
	typedef _Alloc0 _Alloc;
	typedef _Tree_base_types<_Ty, _Alloc> _Myt;

	typedef _Wrap_alloc<_Alloc> _Alty0;
	typedef typename _Alty0::template rebind<_Ty>::other _Alty;

	typedef typename _Get_voidptr<_Alty, typename _Alty::pointer>::type
		_Voidptr;
	typedef _Tree_node<typename _Alty::value_type,
		_Voidptr> _Node;

	typedef typename _Alty::template rebind<_Node>::other _Alnod_type;
	typedef typename _Alnod_type::pointer _Nodeptr;
	typedef _Nodeptr& _Nodepref;

	typedef typename _If<_Is_simple_alloc<_Alty>::value,
		_Tree_simple_types<typename _Alty::value_type>,
		_Tree_iter_types<typename _Alty::value_type,
			typename _Alty::size_type,
			typename _Alty::difference_type,
			typename _Alty::pointer,
			typename _Alty::const_pointer,
			typename _Alty::reference,
			typename _Alty::const_reference,
			_Nodeptr> >::type
		_Val_types;
	};

		
template<class _Val_types>
	class _Tree_val
		: public _Container_base
	{	
public:
	typedef _Tree_val<_Val_types> _Myt;

	typedef typename _Val_types::_Nodeptr _Nodeptr;
	typedef _Nodeptr& _Nodepref;

	typedef typename _Val_types::value_type value_type;
	typedef typename _Val_types::size_type size_type;
	typedef typename _Val_types::difference_type difference_type;
	typedef typename _Val_types::pointer pointer;
	typedef typename _Val_types::const_pointer const_pointer;
	typedef typename _Val_types::reference reference;
	typedef typename _Val_types::const_reference const_reference;

	typedef _Tree_const_iterator<_Myt> const_iterator;
	typedef _Tree_iterator<_Myt> iterator;

	_Tree_val()
		: _Myhead(),
		_Mysize(0)
		{	
		}

	enum _Redbl
		{	
		_Red, _Black};

	static char& _Color(_Nodeptr _Pnode)
		{	
		return ((char&)_Pnode->_Color);
		}

	static char& _Isnil(_Nodeptr _Pnode)
		{	
		return ((char&)_Pnode->_Isnil);
		}

	static _Nodepref _Left(_Nodeptr _Pnode)
		{	
		return ((_Nodepref)_Pnode->_Left);
		}

	static _Nodepref _Parent(_Nodeptr _Pnode)
		{	
		return ((_Nodepref)_Pnode->_Parent);
		}

	static _Nodepref _Right(_Nodeptr _Pnode)
		{	
		return ((_Nodepref)_Pnode->_Right);
		}

	static reference _Myval(_Nodeptr _Pnode)
		{	
		return ((reference)_Pnode->_Myval);
		}

	static _Nodeptr _Max(_Nodeptr _Pnode)
		{	
		while (!_Isnil(_Right(_Pnode)))
			_Pnode = _Right(_Pnode);
		return (_Pnode);
		}

	static _Nodeptr _Min(_Nodeptr _Pnode)
		{	
		while (!_Isnil(_Left(_Pnode)))
			_Pnode = _Left(_Pnode);
		return (_Pnode);
		}

	_Nodeptr _Myhead;	
	size_type _Mysize;	
	};

		
template<class _Traits>
	class _Tree_comp_alloc
	{	
public:
	typedef _Tree_comp_alloc<_Traits> _Myt;

	typedef typename _Traits::allocator_type allocator_type;
	typedef typename _Traits::key_compare key_compare;

	typedef _Tree_base_types<typename _Traits::value_type,
		allocator_type> _Alloc_types;

	typedef typename _Alloc_types::_Alloc _Alloc;
	typedef typename _Alloc_types::_Alnod_type _Alty;
	typedef typename _Alloc_types::_Node _Node;
	typedef typename _Alloc_types::_Nodeptr _Nodeptr;
	typedef typename _Alloc_types::_Val_types _Val_types;

	typedef _Nodeptr& _Nodepref;

	typedef typename _Val_types::value_type value_type;
	typedef typename _Val_types::size_type size_type;
	typedef typename _Val_types::difference_type difference_type;
	typedef typename _Val_types::pointer pointer;
	typedef typename _Val_types::const_pointer const_pointer;
	typedef typename _Val_types::reference reference;
	typedef typename _Val_types::const_reference const_reference;

	typedef _Tree_const_iterator<_Tree_val<_Val_types> > const_iterator;
	typedef _Tree_iterator<_Tree_val<_Val_types> > iterator;

	enum _Redbl
		{	
		_Red, _Black
		};

	static char& _Color(_Nodeptr _Pnode)
		{	
		return (_Tree_val<_Val_types>::_Color(_Pnode));
		}

	static char& _Isnil(_Nodeptr _Pnode)
		{	
		return (_Tree_val<_Val_types>::_Isnil(_Pnode));
		}

	static _Nodepref _Left(_Nodeptr _Pnode)
		{	
		return (_Tree_val<_Val_types>::_Left(_Pnode));
		}

	static _Nodepref _Parent(_Nodeptr _Pnode)
		{	
		return (_Tree_val<_Val_types>::_Parent(_Pnode));
		}

	static _Nodepref _Right(_Nodeptr _Pnode)
		{	
		return (_Tree_val<_Val_types>::_Right(_Pnode));
		}

	static reference _Myval(_Nodeptr _Pnode)
		{	
		return (_Tree_val<_Val_types>::_Myval(_Pnode));
		}

	static _Nodeptr _Max(_Nodeptr _Pnode)
		{	
		return (_Tree_val<_Val_types>::_Max(_Pnode));
		}

	static _Nodeptr _Min(_Nodeptr _Pnode)
		{	
		return (_Tree_val<_Val_types>::_Min(_Pnode));
		}

	_Tree_comp_alloc(const key_compare& _Parg)
		: _Mypair(_One_then_variadic_args_t(), _Parg,
			_Zero_then_variadic_args_t())
		{	
		_Construct();
		}

	template<class _Any_alloc,
		class = enable_if_t<!is_same<decay_t<_Any_alloc>, _Myt>::value> >
		_Tree_comp_alloc(const key_compare& _Parg, _Any_alloc&& _Al)
		: _Mypair(_One_then_variadic_args_t(), _Parg,
			_One_then_variadic_args_t(),
			::std:: forward<_Any_alloc>(_Al))
		{	
		_Construct();
		}

 
	void _Construct()
		{	
		_Myhead() = _Buyheadnode();
		}

	~_Tree_comp_alloc() noexcept
		{	
		_Freeheadnode(_Myhead());
		}

	void _Copy_alloc(const _Alty& _Al)
		{	
		_Pocca(_Getal(), _Al);
		}

	void _Move_alloc(_Alty& _Al)
		{	
		_Pocma(_Getal(), _Al);
		}

 


































































	void _Orphan_all()
		{	
		_Get_data()._Orphan_all();
		}

	void _Swap_all(_Myt& _Right)
		{	
		_Get_data()._Swap_all(_Right._Get_data());
		}

	_Nodeptr _Buyheadnode()
		{	
		_Nodeptr _Pnode = _Getal().allocate(1);

		try {
		_Getal().construct(
			::std:: addressof(_Left(_Pnode)), _Pnode);
		_Getal().construct(
			::std:: addressof(_Parent(_Pnode)), _Pnode);
		_Getal().construct(
			::std:: addressof(_Right(_Pnode)), _Pnode);
		} catch (...) {
		_Getal().deallocate(_Pnode, 1);
		throw;
		}

		_Color(_Pnode) = _Black;
		_Isnil(_Pnode) = true;
		return (_Pnode);
		}

	void _Freeheadnode(_Nodeptr _Pnode)
		{	
		_Getal().destroy(
			::std:: addressof(_Left(_Pnode)));
		_Getal().destroy(
			::std:: addressof(_Parent(_Pnode)));
		_Getal().destroy(
			::std:: addressof(_Right(_Pnode)));
		_Getal().deallocate(_Pnode, 1);
		}

	_Nodeptr _Buynode0()
		{	
		_Nodeptr _Pnode = _Getal().allocate(1);

		try {
		_Getal().construct(
			::std:: addressof(_Left(_Pnode)), _Myhead());
		_Getal().construct(
			::std:: addressof(_Parent(_Pnode)), _Myhead());
		_Getal().construct(
			::std:: addressof(_Right(_Pnode)), _Myhead());
		} catch (...) {
		_Getal().deallocate(_Pnode, 1);
		throw;
		}

		return (_Pnode);
		}

	void _Freenode0(_Nodeptr _Pnode)
		{	
		_Getal().destroy(
			::std:: addressof(_Left(_Pnode)));
		_Getal().destroy(
			::std:: addressof(_Parent(_Pnode)));
		_Getal().destroy(
			::std:: addressof(_Right(_Pnode)));
		_Getal().deallocate(_Pnode, 1);
		}

	template<class... _Valty>
		_Nodeptr _Buynode(_Valty&&... _Val)
		{	
		_Nodeptr _Pnode = _Buynode0();

		this->_Color(_Pnode) = _Red;
		this->_Isnil(_Pnode) = false;

		try {
		this->_Getal().construct(
			::std:: addressof(_Myval(_Pnode)),
				::std:: forward<_Valty>(_Val)...);
		} catch (...) {
		_Freenode0(_Pnode);
		throw;
		}

		return (_Pnode);
		}

	key_compare& _Getcomp() noexcept
		{	
		return (_Mypair._Get_first());
		}

	const key_compare& _Getcomp() const noexcept
		{	
		return (_Mypair._Get_first());
		}

	_Alty& _Getal() noexcept
		{	
		return (_Mypair._Get_second()._Get_first());
		}

	const _Alty& _Getal() const noexcept
		{	
		return (_Mypair._Get_second()._Get_first());
		}

	_Tree_val<_Val_types>& _Get_data() noexcept
		{	
		return (_Mypair._Get_second()._Get_second());
		}

	const _Tree_val<_Val_types>& _Get_data() const noexcept
		{	
		return (_Mypair._Get_second()._Get_second());
		}

	_Nodeptr& _Myhead() noexcept
		{	
		return (_Get_data()._Myhead);
		}

	const _Nodeptr& _Myhead() const noexcept
		{	
		return (_Get_data()._Myhead);
		}

	size_type& _Mysize() noexcept
		{	
		return (_Get_data()._Mysize);
		}

	const size_type& _Mysize() const noexcept
		{	
		return (_Get_data()._Mysize);
		}

private:
	_Compressed_pair<key_compare,
		_Compressed_pair<_Alty, _Tree_val<_Val_types> > > _Mypair;
	};

		
template<class _Traits>
	class _Tree
		: public _Tree_comp_alloc<_Traits>
	{	
public:
	typedef _Tree<_Traits> _Myt;
	typedef _Tree_comp_alloc<_Traits> _Mybase;

	typedef typename _Traits::key_type key_type;
	typedef typename _Traits::value_compare value_compare;
	enum
		{	
		_Multi = _Traits::_Multi};

	typedef typename _Mybase::_Node _Node;
	typedef typename _Mybase::_Nodeptr _Nodeptr;
	typedef typename _Mybase::_Alty _Alty;

	typedef typename _Mybase::key_compare key_compare;
	typedef typename _Mybase::allocator_type allocator_type;

	typedef typename _Mybase::value_type value_type;
	typedef typename _Mybase::size_type size_type;
	typedef typename _Mybase::difference_type difference_type;
	typedef typename _Mybase::pointer pointer;
	typedef typename _Mybase::const_pointer const_pointer;
	typedef typename _Mybase::reference reference;
	typedef typename _Mybase::const_reference const_reference;

	typedef typename _Mybase::const_iterator const_iterator;
	typedef typename _If<is_same<key_type, value_type>::value,
		typename _Mybase::const_iterator,
		typename _Mybase::iterator>::type iterator;

	typedef ::std:: reverse_iterator<iterator> reverse_iterator;
	typedef ::std:: reverse_iterator<const_iterator> const_reverse_iterator;

	typedef pair<iterator, bool> _Pairib;
	typedef pair<iterator, iterator> _Pairii;
	typedef pair<const_iterator, const_iterator> _Paircc;

	struct _Copy_tag
		{	
		};
	struct _Move_tag
		{	
		};

	_Tree(const key_compare& _Parg)
		: _Mybase(_Parg)
		{	
		}

	_Tree(const key_compare& _Parg,
		const allocator_type& _Al)
		: _Mybase(_Parg, _Al)
		{	
		}

	template<class _Any_alloc>
		_Tree(const _Myt& _Right, _Any_alloc&& _Al)
		: _Mybase(_Right.key_comp(), ::std:: forward<_Any_alloc>(_Al))
		{	
		try {
		_Copy(_Right, _Copy_tag());
		} catch (...) {
		_Tidy();
		throw;
		}
		}

	_Tree(_Myt&& _Right)
		: _Mybase(_Right.key_comp(), ::std:: move(_Right._Getal()))
		{	
		_Assign_rv(::std:: forward<_Myt>(_Right), true_type());
		}

	_Tree(_Myt&& _Right, const allocator_type& _Al)
		: _Mybase(_Right.key_comp(), _Al)
		{	
		_Assign_rv(::std:: forward<_Myt>(_Right));
		}

	_Myt& operator=(_Myt&& _Right)
		{	
		if (this != &_Right)
			{	
			clear();
			if (_Alty::propagate_on_container_move_assignment::value
				&& this->_Getal() != _Right._Getal())
				this->_Move_alloc(_Right._Getal());

			_Assign_rv(::std:: forward<_Myt>(_Right));
			}
		return (*this);
		}

	void _Assign_rv(_Myt&& _Right, true_type)
		{	
		this->_Swap_all(_Right);
		_Swap_adl(this->_Getcomp(), _Right._Getcomp());
		_Swap_adl(this->_Myhead(), _Right._Myhead());
		::std:: swap(this->_Mysize(), _Right._Mysize());
		}

	void _Assign_rv(_Myt&& _Right, false_type)
		{	
		if (get_allocator() == _Right.get_allocator())
			_Assign_rv(::std:: forward<_Myt>(_Right), true_type());
		else
			_Copy(_Right, _Move_tag());
		}

	void _Assign_rv(_Myt&& _Right)
		{	
		_Assign_rv(::std:: forward<_Myt>(_Right),
			typename _Alty::propagate_on_container_move_assignment());
		}

	template<class... _Valty>
		_Pairib emplace(_Valty&&... _Val)
		{	
		_Nodeptr _Newnode = this->_Buynode(::std:: forward<_Valty>(_Val)...);
		return (_Insert_nohint(false,
			this->_Myval(_Newnode), _Newnode));
		}

	template<class... _Valty>
		iterator emplace_hint(const_iterator _Where, _Valty&&... _Val)
		{	
		_Nodeptr _Newnode = this->_Buynode(::std:: forward<_Valty>(_Val)...);
		return (_Insert_hint(_Where,
			this->_Myval(_Newnode), _Newnode));
		}

	~_Tree() noexcept
		{	
		_Tidy();
		}

	_Myt& operator=(const _Myt& _Right)
		{	
		if (this != &_Right)
			{	
			clear();
			if (this->_Getal() != _Right._Getal()
				&& _Alty::propagate_on_container_copy_assignment::value)
				this->_Copy_alloc(_Right._Getal());

			this->_Getcomp() = _Right._Getcomp();
			_Copy(_Right, _Copy_tag());
			}
		return (*this);
		}

	iterator begin() noexcept
		{	
		return (iterator(_Lmost(), &this->_Get_data()));
		}

	const_iterator begin() const noexcept
		{	
		return (const_iterator(_Lmost(), &this->_Get_data()));
		}

	iterator end() noexcept
		{	
		return (iterator(this->_Myhead(), &this->_Get_data()));
		}

	const_iterator end() const noexcept
		{	
		return (const_iterator(this->_Myhead(), &this->_Get_data()));
		}

	reverse_iterator rbegin() noexcept
		{	
		return (reverse_iterator(end()));
		}

	const_reverse_iterator rbegin() const noexcept
		{	
		return (const_reverse_iterator(end()));
		}

	reverse_iterator rend() noexcept
		{	
		return (reverse_iterator(begin()));
		}

	const_reverse_iterator rend() const noexcept
		{	
		return (const_reverse_iterator(begin()));
		}

	const_iterator cbegin() const noexcept
		{	
		return (begin());
		}

	const_iterator cend() const noexcept
		{	
		return (end());
		}

	const_reverse_iterator crbegin() const noexcept
		{	
		return (rbegin());
		}

	const_reverse_iterator crend() const noexcept
		{	
		return (rend());
		}

	size_type size() const noexcept
		{	
		return (this->_Mysize());
		}

	size_type max_size() const noexcept
		{	
		return (this->_Getal().max_size());
		}

	bool empty() const noexcept
		{	
		return (size() == 0);
		}

	allocator_type get_allocator() const noexcept
		{	
		allocator_type _Ret(this->_Getal());
		return (_Ret);
		}

	key_compare key_comp() const
		{	
		return (this->_Getcomp());
		}

	value_compare value_comp() const
		{	
		return (value_compare(key_comp()));
		}

	template<bool _Multi2 = _Multi,
		enable_if_t<!_Multi2, int> = 0>
		_Pairib insert(const value_type& _Val)
		{	
		return (_Insert_nohint(false,
			_Val, _Nil()));
		}

	template<bool _Multi2 = _Multi,
		enable_if_t<_Multi2, int> = 0>
		iterator insert(const value_type& _Val)
		{	
		return (_Insert_nohint(false,
			_Val, _Nil()).first);
		}

	template<bool _Multi2 = _Multi,
		enable_if_t<!_Multi2, int> = 0>
		_Pairib insert(value_type&& _Val)
		{	
		return (_Insert_nohint(false,
			::std:: forward<value_type>(_Val), _Nil()));
		}

	template<bool _Multi2 = _Multi,
		enable_if_t<_Multi2, int> = 0>
		iterator insert(value_type&& _Val)
		{	
		return (_Insert_nohint(false,
			::std:: forward<value_type>(_Val), _Nil()).first);
		}

	iterator insert(const_iterator _Where,
		const value_type& _Val)
		{	
		return (_Insert_hint(_Where,
			_Val, _Nil()));
		}

	iterator insert(const_iterator _Where, value_type&& _Val)
		{	
		return (_Insert_hint(_Where,
			::std:: forward<value_type>(_Val), _Nil()));
		}

	template<class _Iter>
		void insert(_Iter _First, _Iter _Last)
		{	
		;
		for (; _First != _Last; ++_First)
			emplace_hint(end(), *_First);
		}

	void insert(::std:: initializer_list<value_type> _Ilist)
		{	
		insert(_Ilist.begin(), _Ilist.end());
		}

	iterator erase(const_iterator _Where)
		{	
 








		_Nodeptr _Erasednode = _Where._Mynode();	
		++_Where;	
 

		_Nodeptr _Fixnode;	
		_Nodeptr _Fixnodeparent;	
		_Nodeptr _Pnode = _Erasednode;

		if (this->_Isnil(this->_Left(_Pnode)))
			_Fixnode = this->_Right(_Pnode);	
		else if (this->_Isnil(this->_Right(_Pnode)))
			_Fixnode = this->_Left(_Pnode);	
		else
			{	
			_Pnode = _Where._Mynode();	
			_Fixnode = this->_Right(_Pnode);	
			}

		if (_Pnode == _Erasednode)
			{	
			_Fixnodeparent = this->_Parent(_Erasednode);
			if (!this->_Isnil(_Fixnode))
				this->_Parent(_Fixnode) = _Fixnodeparent;	

			if (_Root() == _Erasednode)
				_Root() = _Fixnode;	
			else if (this->_Left(_Fixnodeparent) == _Erasednode)
				this->_Left(_Fixnodeparent) = _Fixnode;	
			else
				this->_Right(_Fixnodeparent) =
					_Fixnode;	

			if (_Lmost() == _Erasednode)
				_Lmost() = this->_Isnil(_Fixnode)
					? _Fixnodeparent	
					: this->_Min(_Fixnode);	

			if (_Rmost() == _Erasednode)
				_Rmost() = this->_Isnil(_Fixnode)
					? _Fixnodeparent	
					: this->_Max(_Fixnode);	
			}
		else
			{	
			this->_Parent(this->_Left(_Erasednode)) =
				_Pnode;	
			this->_Left(_Pnode) =
				this->_Left(_Erasednode);	

			if (_Pnode == this->_Right(_Erasednode))
				_Fixnodeparent = _Pnode;	
			else
				{	
				_Fixnodeparent =
					this->_Parent(_Pnode);	
				if (!this->_Isnil(_Fixnode))
					this->_Parent(_Fixnode) = _Fixnodeparent;	
				this->_Left(_Fixnodeparent) = _Fixnode;	
				this->_Right(_Pnode) =
					this->_Right(_Erasednode);	
				this->_Parent(this->_Right(_Erasednode)) =
					_Pnode;	
				}

			if (_Root() == _Erasednode)
				_Root() = _Pnode;	
			else if (this->_Left(this->_Parent(_Erasednode)) == _Erasednode)
				this->_Left(this->_Parent(_Erasednode)) =
					_Pnode;	
			else
				this->_Right(this->_Parent(_Erasednode)) =
					_Pnode;	

			this->_Parent(_Pnode) =
				this->_Parent(_Erasednode);	
			::std:: swap(this->_Color(_Pnode),
				this->_Color(_Erasednode));	
			}

		if (this->_Color(_Erasednode) == this->_Black)
			{	
			for (; _Fixnode != _Root()
				&& this->_Color(_Fixnode) == this->_Black;
				_Fixnodeparent = this->_Parent(_Fixnode))
				if (_Fixnode == this->_Left(_Fixnodeparent))
					{	
					_Pnode = this->_Right(_Fixnodeparent);
					if (this->_Color(_Pnode) == this->_Red)
						{	
						this->_Color(_Pnode) = this->_Black;
						this->_Color(_Fixnodeparent) = this->_Red;
						_Lrotate(_Fixnodeparent);
						_Pnode = this->_Right(_Fixnodeparent);
						}

					if (this->_Isnil(_Pnode))
						_Fixnode = _Fixnodeparent;	
					else if (this->_Color(this->_Left(_Pnode)) == this->_Black
						&& this->_Color(this->_Right(_Pnode)) == this->_Black)
						{	
						this->_Color(_Pnode) = this->_Red;
						_Fixnode = _Fixnodeparent;
						}
					else
						{	
						if (this->_Color(this->_Right(_Pnode))
							== this->_Black)
							{	
							this->_Color(this->_Left(_Pnode)) = this->_Black;
							this->_Color(_Pnode) = this->_Red;
							_Rrotate(_Pnode);
							_Pnode = this->_Right(_Fixnodeparent);
							}

						this->_Color(_Pnode) = this->_Color(_Fixnodeparent);
						this->_Color(_Fixnodeparent) = this->_Black;
						this->_Color(this->_Right(_Pnode)) = this->_Black;
						_Lrotate(_Fixnodeparent);
						break;	
						}
					}
				else
					{	
					_Pnode = this->_Left(_Fixnodeparent);
					if (this->_Color(_Pnode) == this->_Red)
						{	
						this->_Color(_Pnode) = this->_Black;
						this->_Color(_Fixnodeparent) = this->_Red;
						_Rrotate(_Fixnodeparent);
						_Pnode = this->_Left(_Fixnodeparent);
						}

					if (this->_Isnil(_Pnode))
						_Fixnode = _Fixnodeparent;	
					else if (this->_Color(this->_Right(_Pnode)) ==
						this->_Black
						&& this->_Color(this->_Left(_Pnode)) == this->_Black)
						{	
						this->_Color(_Pnode) = this->_Red;
						_Fixnode = _Fixnodeparent;
						}
					else
						{	
						if (this->_Color(this->_Left(_Pnode)) == this->_Black)
							{	
							this->_Color(this->_Right(_Pnode)) = this->_Black;
							this->_Color(_Pnode) = this->_Red;
							_Lrotate(_Pnode);
							_Pnode = this->_Left(_Fixnodeparent);
							}

						this->_Color(_Pnode) = this->_Color(_Fixnodeparent);
						this->_Color(_Fixnodeparent) = this->_Black;
						this->_Color(this->_Left(_Pnode)) = this->_Black;
						_Rrotate(_Fixnodeparent);
						break;	
						}
					}

			this->_Color(_Fixnode) = this->_Black;	
			}

		this->_Getal().destroy(
			::std:: addressof(this->_Myval(_Erasednode)));	

		this->_Getal().deallocate(_Erasednode, 1);

		if (0 < this->_Mysize())
			--this->_Mysize();

		return (iterator(_Where._Ptr,
			&this->_Get_data()));	
		}

	iterator erase(const_iterator _First, const_iterator _Last)
		{	
		if (_First == begin() && _Last == end())
			{	
			clear();
			return (begin());
			}
		else
			{	
			while (_First != _Last)
				erase(_First++);
			return (iterator(_First._Ptr, &this->_Get_data()));
			}
		}

	size_type erase(const key_type& _Keyval)
		{	
		_Pairii _Where = equal_range(_Keyval);
		size_type _Num = ::std:: distance(_Where.first, _Where.second);
		erase(_Where.first, _Where.second);
		return (_Num);
		}

	void clear() noexcept
		{	
 



		_Erase(_Root());
		_Root() = this->_Myhead();
		_Lmost() = this->_Myhead();
		_Rmost() = this->_Myhead();
		this->_Mysize() = 0;
		}

	iterator find(const key_type& _Keyval)
		{	
		iterator _Where = lower_bound(_Keyval);
		return (_Where == end()
			|| this->_Getcomp()(_Keyval, this->_Key(_Where._Mynode()))
					? end() : _Where);
		}

	const_iterator find(const key_type& _Keyval) const
		{	
		const_iterator _Where = lower_bound(_Keyval);
		return (_Where == end()
			|| this->_Getcomp()(_Keyval, this->_Key(_Where._Mynode()))
					? end() : _Where);
		}

	template<class _Other,
		class _Mycomp = key_compare,
		class = typename _Mycomp::is_transparent>
		iterator find(const _Other& _Keyval)
		{	
		iterator _Where = lower_bound(_Keyval);
		return (_Where == end()
			|| this->_Getcomp()(_Keyval, this->_Key(_Where._Mynode()))
					? end() : _Where);
		}

	template<class _Other,
		class _Mycomp = key_compare,
		class = typename _Mycomp::is_transparent>
		const_iterator find(const _Other& _Keyval) const
		{	
		const_iterator _Where = lower_bound(_Keyval);
		return (_Where == end()
			|| this->_Getcomp()(_Keyval, this->_Key(_Where._Mynode()))
					? end() : _Where);
		}

	size_type count(const key_type& _Keyval) const
		{	
		_Paircc _Ans = equal_range(_Keyval);
		return (::std:: distance(_Ans.first, _Ans.second));
		}

	template<class _Other,
		class _Mycomp = key_compare,
		class = typename _Mycomp::is_transparent>
		size_type count(const _Other& _Keyval) const
		{	
		_Paircc _Ans = equal_range(_Keyval);
		return (::std:: distance(_Ans.first, _Ans.second));
		}

	iterator lower_bound(const key_type& _Keyval)
		{	
		return (iterator(_Lbound(_Keyval), &this->_Get_data()));
		}

	const_iterator lower_bound(const key_type& _Keyval) const
		{	
		return (const_iterator(_Lbound(_Keyval), &this->_Get_data()));
		}

	template<class _Other,
		class _Mycomp = key_compare,
		class = typename _Mycomp::is_transparent>
		iterator lower_bound(const _Other& _Keyval)
		{	
		return (iterator(_Lbound(_Keyval), &this->_Get_data()));
		}

	template<class _Other,
		class _Mycomp = key_compare,
		class = typename _Mycomp::is_transparent>
		const_iterator lower_bound(const _Other& _Keyval) const
		{	
		return (const_iterator(_Lbound(_Keyval), &this->_Get_data()));
		}

	iterator upper_bound(const key_type& _Keyval)
		{	
		return (iterator(_Ubound(_Keyval), &this->_Get_data()));
		}

	const_iterator upper_bound(const key_type& _Keyval) const
		{	
		return (const_iterator(_Ubound(_Keyval), &this->_Get_data()));
		}

	template<class _Other,
		class _Mycomp = key_compare,
		class = typename _Mycomp::is_transparent>
		iterator upper_bound(const _Other& _Keyval)
		{	
		return (iterator(_Ubound(_Keyval), &this->_Get_data()));
		}

	template<class _Other,
		class _Mycomp = key_compare,
		class = typename _Mycomp::is_transparent>
		const_iterator upper_bound(const _Other& _Keyval) const
		{	
		return (const_iterator(_Ubound(_Keyval), &this->_Get_data()));
		}

	_Pairii equal_range(const key_type& _Keyval)
		{	
		return (_Eqrange(_Keyval));
		}

	_Paircc equal_range(const key_type& _Keyval) const
		{	
		return (_Eqrange(_Keyval));
		}

	template<class _Other,
		class _Mycomp = key_compare,
		class = typename _Mycomp::is_transparent>
		_Pairii equal_range(const _Other& _Keyval)
		{	
		return (_Eqrange(_Keyval));
		}

	template<class _Other,
		class _Mycomp = key_compare,
		class = typename _Mycomp::is_transparent>
		_Paircc equal_range(const _Other& _Keyval) const
		{	
		return (_Eqrange(_Keyval));
		}

	void swap(_Myt& _Right)
		{	
		if (this != &_Right)
			{	
			_Pocs(this->_Getal(), _Right._Getal());
			this->_Swap_all(_Right);
			_Swap_adl(this->_Getcomp(), _Right._Getcomp());
			_Swap_adl(this->_Myhead(), _Right._Myhead());
			::std:: swap(this->_Mysize(), _Right._Mysize());
			}
		}

protected:
	template<class _Valty>
		_Nodeptr _Buynode_if_nil(_Nodeptr _Node, _Valty&&)
		{	
		return (_Node);
		}

	template<class _Valty>
		_Nodeptr _Buynode_if_nil(_Nil, _Valty&& _Val)
		{	
		return (this->_Buynode(::std:: forward<_Valty>(_Val)));
		}

	void _Destroy_if_not_nil(_Nodeptr _Newnode)
		{	
		this->_Getal().destroy(
			::std:: addressof(this->_Myval(_Newnode)));

		this->_Getal().deallocate(_Newnode, 1);
		}

	void _Destroy_if_not_nil(_Nil)
		{	
		}

	template<class _Valty,
		class _Nodety>
		iterator _Insert_hint(const_iterator _Where,
			_Valty&& _Val, _Nodety _Newnode)
		{	
		const_iterator _Next;
		bool _Leftish = false;	

		try {

 




		if (size() == 0)
			return (_Insert_at(true, this->_Myhead(),
				::std:: forward<_Valty>(_Val), _Newnode));	
		else if (this->_Multi)
			{	
			if (_Where == begin())
				{	
				if (!this->_Getcomp()(this->_Key(_Where._Mynode()), this->_Kfn(_Val)))
					return (_Insert_at(true, _Where._Mynode(),
						::std:: forward<_Valty>(_Val), _Newnode));
				_Leftish = true;	
				}
			else if (_Where == end())
				{	
				if (!this->_Getcomp()(this->_Kfn(_Val), this->_Key(_Rmost())))
					return (_Insert_at(false, _Rmost(),
						::std:: forward<_Valty>(_Val), _Newnode));
				}
			else if (!this->_Getcomp()(this->_Key(_Where._Mynode()), this->_Kfn(_Val))
				&& !this->_Getcomp()(this->_Kfn(_Val), this->_Key((--(_Next = _Where))._Mynode())))
				{	
				if (this->_Isnil(this->_Right(_Next._Mynode())))
					return (_Insert_at(false, _Next._Mynode(),
						::std:: forward<_Valty>(_Val), _Newnode));
				else
					return (_Insert_at(true, _Where._Mynode(),
						::std:: forward<_Valty>(_Val), _Newnode));
				}
			else if (!this->_Getcomp()(this->_Kfn(_Val), this->_Key(_Where._Mynode()))
				&& (++(_Next = _Where) == end()
					|| !this->_Getcomp()(this->_Key(_Next._Mynode()), this->_Kfn(_Val))))
				{	
				if (this->_Isnil(this->_Right(_Where._Mynode())))
					return (_Insert_at(false, _Where._Mynode(),
						::std:: forward<_Valty>(_Val), _Newnode));
				else
					return (_Insert_at(true, _Next._Mynode(),
						::std:: forward<_Valty>(_Val), _Newnode));
				}
			else
				_Leftish = true;	
			}
		else
			{	
			if (_Where == begin())
				{	
				if (this->_Getcomp()(this->_Kfn(_Val), this->_Key(_Where._Mynode())))
					return (_Insert_at(true, _Where._Mynode(),
						::std:: forward<_Valty>(_Val), _Newnode));
				}
			else if (_Where == end())
				{	
				if (this->_Getcomp()(this->_Key(_Rmost()), this->_Kfn(_Val)))
					return (_Insert_at(false, _Rmost(),
						::std:: forward<_Valty>(_Val), _Newnode));
				}
			else if (this->_Getcomp()(this->_Kfn(_Val), this->_Key(_Where._Mynode()))
				&& this->_Getcomp()(this->_Key((--(_Next = _Where))._Mynode()), this->_Kfn(_Val)))
				{	
				if (this->_Isnil(this->_Right(_Next._Mynode())))
					return (_Insert_at(false, _Next._Mynode(),
						::std:: forward<_Valty>(_Val), _Newnode));
				else
					return (_Insert_at(true, _Where._Mynode(),
						::std:: forward<_Valty>(_Val), _Newnode));
				}
			else if (this->_Getcomp()(this->_Key(_Where._Mynode()), this->_Kfn(_Val))
				&& (++(_Next = _Where) == end()
					|| this->_Getcomp()(this->_Kfn(_Val), this->_Key(_Next._Mynode()))))
				{	
				if (this->_Isnil(this->_Right(_Where._Mynode())))
					return (_Insert_at(false, _Where._Mynode(),
						::std:: forward<_Valty>(_Val), _Newnode));
				else
					return (_Insert_at(true, _Next._Mynode(),
						::std:: forward<_Valty>(_Val), _Newnode));
				}
			}
		} catch (...) {
		_Destroy_if_not_nil(_Newnode);
		throw;
		}

		return (_Insert_nohint(_Leftish,
			::std:: forward<_Valty>(_Val), _Newnode).first);
		}

	template<class _Valty,
		class _Nodety>
		_Pairib _Insert_nohint(bool _Leftish,
			_Valty&& _Val, _Nodety _Newnode)
		{	
		try {
		_Nodeptr _Trynode = _Root();
		_Nodeptr _Wherenode = this->_Myhead();
		bool _Addleft = true;	

		while (!this->_Isnil(_Trynode))
			{	
			_Wherenode = _Trynode;
			if (_Leftish)
				_Addleft = !this->_Getcomp()(this->_Key(_Trynode), this->_Kfn(_Val));	
			else
				_Addleft = this->_Getcomp()(this->_Kfn(_Val), this->_Key(_Trynode));	
			_Trynode = _Addleft ? this->_Left(_Trynode)
				: this->_Right(_Trynode);
			}

		if (this->_Multi)
			return (_Pairib(_Insert_at(_Addleft, _Wherenode,
				::std:: forward<_Valty>(_Val), _Newnode), true));
		else
			{	
			iterator _Where = iterator(_Wherenode, &this->_Get_data());
			if (!_Addleft)
				;	
			else if (_Where == begin())
				return (_Pairib(_Insert_at(true, _Wherenode,
					::std:: forward<_Valty>(_Val), _Newnode), true));
			else
				--_Where;	

			if (this->_Getcomp()(this->_Key(_Where._Mynode()), this->_Kfn(_Val)))
				return (_Pairib(_Insert_at(_Addleft, _Wherenode,
					::std:: forward<_Valty>(_Val), _Newnode), true));
			else
				{	
				_Destroy_if_not_nil(_Newnode);
				return (_Pairib(_Where, false));
				}
			}
		} catch (...) {
		_Destroy_if_not_nil(_Newnode);
		throw;
		}
		}

	template<class _Valty,
		class _Nodety>
		iterator _Insert_at(bool _Addleft, _Nodeptr _Wherenode,
		_Valty&& _Val, _Nodety _Node)
		{	
		if (max_size() - 1 <= this->_Mysize())
			{	
			_Destroy_if_not_nil(_Node);
			_Xlength_error("map/set<T> too long");
			}
		_Nodeptr _Newnode = _Buynode_if_nil(_Node,
			::std:: forward<_Valty>(_Val));

		++this->_Mysize();
		_Newnode->_Parent = _Wherenode;

		if (_Wherenode == this->_Myhead())
			{	
			_Root() = _Newnode;
			_Lmost() = _Newnode;
			_Rmost() = _Newnode;
			}
		else if (_Addleft)
			{	
			this->_Left(_Wherenode) = _Newnode;
			if (_Wherenode == _Lmost())
				_Lmost() = _Newnode;
			}
		else
			{	
			this->_Right(_Wherenode) = _Newnode;
			if (_Wherenode == _Rmost())
				_Rmost() = _Newnode;
			}

		for (_Nodeptr _Pnode = _Newnode;
			this->_Color(this->_Parent(_Pnode)) == this->_Red; )
			if (this->_Parent(_Pnode)
				== this->_Left(this->_Parent(this->_Parent(_Pnode))))
				{	
				_Wherenode =
					this->_Right(this->_Parent(this->_Parent(_Pnode)));
				if (this->_Color(_Wherenode) == this->_Red)
					{	
					this->_Color(this->_Parent(_Pnode)) = this->_Black;
					this->_Color(_Wherenode) = this->_Black;
					this->_Color(this->_Parent(this->_Parent(_Pnode)))
						= this->_Red;
					_Pnode = this->_Parent(this->_Parent(_Pnode));
					}
				else
					{	
					if (_Pnode == this->_Right(this->_Parent(_Pnode)))
						{	
						_Pnode = this->_Parent(_Pnode);
						_Lrotate(_Pnode);
						}
					this->_Color(this->_Parent(_Pnode)) =
						this->_Black;	
					this->_Color(this->_Parent(this->_Parent(_Pnode))) =
						this->_Red;
					_Rrotate(this->_Parent(this->_Parent(_Pnode)));
					}
				}
			else
				{	
				_Wherenode =
					this->_Left(this->_Parent(this->_Parent(_Pnode)));
				if (this->_Color(_Wherenode) == this->_Red)
					{	
					this->_Color(this->_Parent(_Pnode)) = this->_Black;
					this->_Color(_Wherenode) = this->_Black;
					this->_Color(this->_Parent(this->_Parent(_Pnode))) =
						this->_Red;
					_Pnode = this->_Parent(this->_Parent(_Pnode));
					}
				else
					{	
					if (_Pnode == this->_Left(this->_Parent(_Pnode)))
						{	
						_Pnode = this->_Parent(_Pnode);
						_Rrotate(_Pnode);
						}
					this->_Color(this->_Parent(_Pnode)) =
						this->_Black;	
					this->_Color(this->_Parent(this->_Parent(_Pnode))) =
						this->_Red;
					_Lrotate(this->_Parent(this->_Parent(_Pnode)));
					}
				}

		this->_Color(_Root()) = this->_Black;	
		return (iterator(_Newnode, &this->_Get_data()));
		}

	template<class _Moveit>
		void _Copy(const _Myt& _Right, _Moveit _Movefl)
		{	
		_Root() = _Copy_nodes(_Right._Root(), this->_Myhead(), _Movefl);
		this->_Mysize() = _Right.size();
		if (!this->_Isnil(_Root()))
			{	
			_Lmost() = this->_Min(_Root());
			_Rmost() = this->_Max(_Root());
			}
		else
			{	
			_Lmost() = this->_Myhead();
			_Rmost() = this->_Myhead();
			}
		}

	template<class _Ty,
		class _Is_set>
		_Nodeptr _Copy_or_move(_Ty& _Val, _Copy_tag, _Is_set)
		{	
		return (this->_Buynode(_Val));
		}

	template<class _Ty>
		_Nodeptr _Copy_or_move(_Ty& _Val, _Move_tag, true_type)
		{	
		return (this->_Buynode(::std:: move(_Val)));
		}

	template<class _Ty>
		_Nodeptr _Copy_or_move(_Ty& _Val, _Move_tag, false_type)
		{	
		return (this->_Buynode(
			::std:: move(const_cast<key_type&>(_Val.first)),
			::std:: move(_Val.second)));
		}

	template<class _Moveit>
		_Nodeptr _Copy_nodes(_Nodeptr _Rootnode, _Nodeptr _Wherenode,
			_Moveit _Movefl)
		{	
		_Nodeptr _Newroot = this->_Myhead();	

		if (!this->_Isnil(_Rootnode))
			{	
			typename is_same<key_type, value_type>::type _Is_set;
			_Nodeptr _Pnode = _Copy_or_move(
				this->_Myval(_Rootnode), _Movefl, _Is_set);
			_Pnode->_Parent = _Wherenode;
			_Pnode->_Color = this->_Color(_Rootnode);
			if (this->_Isnil(_Newroot))
				_Newroot = _Pnode;	

			try {
			this->_Left(_Pnode) =
				_Copy_nodes(this->_Left(_Rootnode), _Pnode, _Movefl);
			this->_Right(_Pnode) =
				_Copy_nodes(this->_Right(_Rootnode), _Pnode, _Movefl);
			} catch (...) {
			_Erase(_Newroot);	
			throw;
			}
			}

		return (_Newroot);	
		}

	template<class _Other>
		_Paircc _Eqrange(const _Other& _Keyval) const
		{	
		_Nodeptr _Pnode = _Root();
		_Nodeptr _Lonode = this->_Myhead();	
		_Nodeptr _Hinode = this->_Myhead();	

		while (!this->_Isnil(_Pnode))
			if (this->_Getcomp()(this->_Key(_Pnode), _Keyval))
				_Pnode = this->_Right(_Pnode);	
			else
				{	
				if (this->_Isnil(_Hinode)
						&& this->_Getcomp()(_Keyval, this->_Key(_Pnode)))
					_Hinode = _Pnode;	
				_Lonode = _Pnode;
				_Pnode = this->_Left(_Pnode);	
				}

		_Pnode = this->_Isnil(_Hinode) ? _Root()
			: this->_Left(_Hinode);	
		while (!this->_Isnil(_Pnode))
			if (this->_Getcomp()(_Keyval, this->_Key(_Pnode)))
				{	
				_Hinode = _Pnode;
				_Pnode = this->_Left(_Pnode);	
				}
			else
				_Pnode = this->_Right(_Pnode);	

		const_iterator _First = const_iterator(_Lonode, &this->_Get_data());
		const_iterator _Last = const_iterator(_Hinode, &this->_Get_data());
		return (_Paircc(_First, _Last));
		}

	template<class _Other>
		_Pairii _Eqrange(const _Other& _Keyval)
		{	
		_Paircc _Ans(static_cast<const _Myt *>(this)->_Eqrange(_Keyval));
		iterator _First = iterator(_Ans.first._Ptr, &this->_Get_data());
		iterator _Last = iterator(_Ans.second._Ptr, &this->_Get_data());
		return (_Pairii(_First, _Last));
		}

	void _Erase(_Nodeptr _Rootnode)
		{	
		for (_Nodeptr _Pnode = _Rootnode;
			!this->_Isnil(_Pnode); _Rootnode = _Pnode)
			{	
			_Erase(this->_Right(_Pnode));
			_Pnode = this->_Left(_Pnode);
			this->_Getal().destroy(
				::std:: addressof(this->_Myval(_Rootnode)));

			this->_Getal().deallocate(_Rootnode, 1);
			}
		}

	bool _Compare(const key_type& _Left, const key_type& _Right) const
		{	
		return (this->_Getcomp()(_Left, _Right));
		}

	template<class _Ty1,
		class _Ty2>
		bool _Compare(const _Ty1& _Left, const _Ty2& _Right) const
		{	
		return (this->_Getcomp()(_Left, _Right));
		}

	template<class _Other>
		_Nodeptr _Lbound(const _Other& _Keyval) const
		{	
		_Nodeptr _Pnode = _Root();
		_Nodeptr _Wherenode = this->_Myhead();	

		while (!this->_Isnil(_Pnode))
			if (_Compare(this->_Key(_Pnode), _Keyval))
				_Pnode = this->_Right(_Pnode);	
			else
				{	
				_Wherenode = _Pnode;
				_Pnode = this->_Left(_Pnode);	
				}

		return (_Wherenode);	
		}

	_Nodeptr& _Lmost() const
		{	
		return (this->_Left(this->_Myhead()));
		}

	void _Lrotate(_Nodeptr _Wherenode)
		{	
		_Nodeptr _Pnode = this->_Right(_Wherenode);
		this->_Right(_Wherenode) = this->_Left(_Pnode);

		if (!this->_Isnil(this->_Left(_Pnode)))
			this->_Parent(this->_Left(_Pnode)) = _Wherenode;
		this->_Parent(_Pnode) = this->_Parent(_Wherenode);

		if (_Wherenode == _Root())
			_Root() = _Pnode;
		else if (_Wherenode == this->_Left(this->_Parent(_Wherenode)))
			this->_Left(this->_Parent(_Wherenode)) = _Pnode;
		else
			this->_Right(this->_Parent(_Wherenode)) = _Pnode;

		this->_Left(_Pnode) = _Wherenode;
		this->_Parent(_Wherenode) = _Pnode;
		}

	_Nodeptr& _Rmost() const
		{	
		return (this->_Right(this->_Myhead()));
		}

	_Nodeptr& _Root() const
		{	
		return (this->_Parent(this->_Myhead()));
		}

	void _Rrotate(_Nodeptr _Wherenode)
		{	
		_Nodeptr _Pnode = this->_Left(_Wherenode);
		this->_Left(_Wherenode) = this->_Right(_Pnode);

		if (!this->_Isnil(this->_Right(_Pnode)))
			this->_Parent(this->_Right(_Pnode)) = _Wherenode;
		this->_Parent(_Pnode) = this->_Parent(_Wherenode);

		if (_Wherenode == _Root())
			_Root() = _Pnode;
		else if (_Wherenode == this->_Right(this->_Parent(_Wherenode)))
			this->_Right(this->_Parent(_Wherenode)) = _Pnode;
		else
			this->_Left(this->_Parent(_Wherenode)) = _Pnode;

		this->_Right(_Pnode) = _Wherenode;
		this->_Parent(_Wherenode) = _Pnode;
		}

	template<class _Other>
		_Nodeptr _Ubound(const _Other& _Keyval) const
		{	
		_Nodeptr _Pnode = _Root();
		_Nodeptr _Wherenode = this->_Myhead();	

		while (!this->_Isnil(_Pnode))
			if (_Compare(_Keyval, this->_Key(_Pnode)))
				{	
				_Wherenode = _Pnode;
				_Pnode = this->_Left(_Pnode);	
				}
			else
				_Pnode = this->_Right(_Pnode);	

		return (_Wherenode);	
		}

 

















	void _Tidy()
		{	
		erase(begin(), end());
		}

	const key_type& _Kfn(const value_type& _Val) const
		{	
		return (_Traits::_Kfn(_Val));
		}

	const key_type& _Key(_Nodeptr _Pnode) const
		{	
		return ((const key_type&)this->_Kfn(this->_Myval(_Pnode)));
		}
	};

		
template<class _Traits> inline
	bool operator==(const _Tree<_Traits>& _Left, const _Tree<_Traits>& _Right)
	{	
	return (_Left.size() == _Right.size()
		&& ::std:: equal(_Left.begin(), _Left.end(), _Right.begin()));
	}

template<class _Traits> inline
	bool operator!=(const _Tree<_Traits>& _Left, const _Tree<_Traits>& _Right)
	{	
	return (!(_Left == _Right));
	}

template<class _Traits> inline
	bool operator<(const _Tree<_Traits>& _Left, const _Tree<_Traits>& _Right)
	{	
	return (::std:: lexicographical_compare(_Left.begin(), _Left.end(),
		_Right.begin(), _Right.end()));
	}

template<class _Traits> inline
	bool operator>(const _Tree<_Traits>& _Left, const _Tree<_Traits>& _Right)
	{	
	return (_Right < _Left);
	}

template<class _Traits> inline
	bool operator<=(const _Tree<_Traits>& _Left, const _Tree<_Traits>& _Right)
	{	
	return (!(_Right < _Left));
	}

template<class _Traits> inline
	bool operator>=(const _Tree<_Traits>& _Left, const _Tree<_Traits>& _Right)
	{	
	return (!(_Left < _Right));
	}
}

 
 #pragma warning(pop)
 #pragma pack(pop)










 #pragma pack(push,8)
 #pragma warning(push,3)
 
 
namespace std {
		
template<class _Kty,	
	class _Ty,	
	class _Pr,	
	class _Alloc,	
	bool _Mfl>	
	class _Tmap_traits
	{	
public:
	typedef _Kty key_type;
	typedef pair<const _Kty, _Ty> value_type;
	typedef _Pr key_compare;
	typedef _Alloc allocator_type;

	enum
		{	
		_Multi = _Mfl};

	class value_compare
		{	
		friend class _Tmap_traits<_Kty, _Ty, _Pr, _Alloc, _Mfl>;

	public:
		typedef value_type first_argument_type;
		typedef value_type second_argument_type;
		typedef bool result_type;

		bool operator()(const value_type& _Left,
			const value_type& _Right) const
			{	
			return (comp(_Left.first, _Right.first));
			}

		value_compare(key_compare _Pred)
			: comp(_Pred)
			{	
			}

	protected:
		key_compare comp;	
		};

	template<class _Ty1,
		class _Ty2>
		static const _Kty& _Kfn(const pair<_Ty1, _Ty2>& _Val)
		{	
		return (_Val.first);
		}
	};

		
template<class _Kty,
	class _Ty,
	class _Pr = less<_Kty>,
	class _Alloc = allocator<pair<const _Kty, _Ty> > >
	class map
		: public _Tree<_Tmap_traits<_Kty, _Ty, _Pr, _Alloc, false> >
	{	
public:
	typedef map<_Kty, _Ty, _Pr, _Alloc> _Myt;
	typedef _Tree<_Tmap_traits<_Kty, _Ty, _Pr, _Alloc, false> > _Mybase;
	typedef _Kty key_type;
	typedef _Ty mapped_type;
	typedef _Pr key_compare;
	typedef typename _Mybase::value_compare value_compare;
	typedef typename _Mybase::allocator_type allocator_type;
	typedef typename _Mybase::size_type size_type;
	typedef typename _Mybase::difference_type difference_type;
	typedef typename _Mybase::pointer pointer;
	typedef typename _Mybase::const_pointer const_pointer;
	typedef typename _Mybase::reference reference;
	typedef typename _Mybase::const_reference const_reference;
	typedef typename _Mybase::iterator iterator;
	typedef typename _Mybase::const_iterator const_iterator;
	typedef typename _Mybase::reverse_iterator reverse_iterator;
	typedef typename _Mybase::const_reverse_iterator
		const_reverse_iterator;
	typedef typename _Mybase::value_type value_type;

	typedef typename _Mybase::_Alty _Alty;
	typedef typename _Mybase::_Pairib _Pairib;

	map()
		: _Mybase(key_compare())
		{	
		}

	explicit map(const allocator_type& _Al)
		: _Mybase(key_compare(), _Al)
		{	
		}

	map(const _Myt& _Right)
		: _Mybase(_Right,
			_Right._Getal().select_on_container_copy_construction())
		{	
		}

	map(const _Myt& _Right, const allocator_type& _Al)
		: _Mybase(_Right, _Al)
		{	
		}

	explicit map(const key_compare& _Pred)
		: _Mybase(_Pred)
		{	
		}

	map(const key_compare& _Pred, const allocator_type& _Al)
		: _Mybase(_Pred, _Al)
		{	
		}

	template<class _Iter>
		map(_Iter _First, _Iter _Last)
		: _Mybase(key_compare())
		{	
		insert(_First, _Last);
		}

	template<class _Iter>
		map(_Iter _First, _Iter _Last,
			const key_compare& _Pred)
		: _Mybase(_Pred)
		{	
		insert(_First, _Last);
		}

	template<class _Iter>
		map(_Iter _First, _Iter _Last,
			const key_compare& _Pred, const allocator_type& _Al)
		: _Mybase(_Pred, _Al)
		{	
		insert(_First, _Last);
		}

	_Myt& operator=(const _Myt& _Right)
		{	
		_Mybase::operator=(_Right);
		return (*this);
		}

	map(_Myt&& _Right)
		: _Mybase(::std:: move(_Right))
		{	
		}

	map(_Myt&& _Right, const allocator_type& _Al)
		: _Mybase(::std:: move(_Right), _Al)
		{	
		}

	_Myt& operator=(_Myt&& _Right)
		noexcept(_Alty::is_always_equal::value && is_nothrow_move_assignable<_Pr>::value)
		{	
		_Mybase::operator=(::std:: move(_Right));
		return (*this);
		}

	mapped_type& operator[](key_type&& _Keyval)
		{	
		return (try_emplace(::std:: move(_Keyval)).first->second);
		}

	void swap(_Myt& _Right)
		noexcept(_Alty::is_always_equal::value && _Is_nothrow_swappable<_Pr>::value)
		{	
		_Mybase::swap(_Right);
		}

	using _Mybase::insert;

	template<class _Valty,
		class = enable_if_t<is_constructible<value_type, _Valty>::value> >
		_Pairib insert(_Valty&& _Val)
		{	
		return (this->emplace(::std:: forward<_Valty>(_Val)));
		}

	template<class _Valty,
		class = enable_if_t<is_constructible<value_type, _Valty>::value> >
		iterator insert(const_iterator _Where, _Valty&& _Val)
		{	
		return (this->emplace_hint(_Where, ::std:: forward<_Valty>(_Val)));
		}

	template<class _Keyty,
		class... _Mappedty>
		_Pairib _Try_emplace(_Keyty&& _Keyval,
			_Mappedty&&... _Mapval)
		{	
		iterator _Where = _Mybase::lower_bound(_Keyval);
		if (_Where == _Mybase::end()
			|| _Mybase::_Getcomp()(_Keyval, _Mybase::_Key(_Where._Mynode())))
			return (_Pairib(
				_Mybase::emplace_hint(_Where,
					piecewise_construct,
					::std:: forward_as_tuple(
						::std:: forward<_Keyty>(_Keyval)),
					::std:: forward_as_tuple(
						::std:: forward<_Mappedty>(_Mapval)...)),
				true));
		else
			return (_Pairib(_Where, false));
		}

	template<class... _Mappedty>
		_Pairib try_emplace(const key_type& _Keyval,
			_Mappedty&&... _Mapval)
		{	
		return (_Try_emplace(_Keyval, ::std:: forward<_Mappedty>(_Mapval)...));
		}

	template<class... _Mappedty>
		iterator try_emplace(const_iterator, const key_type& _Keyval,
			_Mappedty&&... _Mapval)
		{	
		return (_Try_emplace(_Keyval,
			::std:: forward<_Mappedty>(_Mapval)...).first);
		}

	template<class... _Mappedty>
		_Pairib try_emplace(key_type&& _Keyval,
			_Mappedty&&... _Mapval)
		{	
		return (_Try_emplace(::std:: move(_Keyval),
			::std:: forward<_Mappedty>(_Mapval)...));
		}

	template<class... _Mappedty>
		iterator try_emplace(const_iterator, key_type&& _Keyval,
			_Mappedty&&... _Mapval)
		{	
		return (_Try_emplace(::std:: move(_Keyval),
			::std:: forward<_Mappedty>(_Mapval)...).first);
		}

	template<class _Keyty,
		class _Mappedty>
		_Pairib _Insert_or_assign(_Keyty&& _Keyval,
			_Mappedty&& _Mapval)
		{	
		iterator _Where = _Mybase::lower_bound(_Keyval);
		if (_Where == _Mybase::end()
			|| _Mybase::_Getcomp()(_Keyval, _Mybase::_Key(_Where._Mynode())))
			return (_Pairib(
				_Mybase::emplace_hint(_Where,
					::std:: forward<_Keyty>(_Keyval),
					::std:: forward<_Mappedty>(_Mapval)),
				true));
		else
			{	
			_Where->second = ::std:: forward<_Mappedty>(_Mapval);
			return (_Pairib(_Where, false));
			}
		}

	template<class _Mappedty>
		_Pairib insert_or_assign(const key_type& _Keyval,
			_Mappedty&& _Mapval)
		{	
		return (_Insert_or_assign(_Keyval,
			::std:: forward<_Mappedty>(_Mapval)));
		}

	template<class _Mappedty>
		iterator insert_or_assign(const_iterator, const key_type& _Keyval,
			_Mappedty&& _Mapval)
		{	
		return (_Insert_or_assign(_Keyval,
			::std:: forward<_Mappedty>(_Mapval)).first);
		}

	template<class _Mappedty>
		_Pairib insert_or_assign(key_type&& _Keyval,
			_Mappedty&& _Mapval)
		{	
		return (_Insert_or_assign(::std:: move(_Keyval),
			::std:: forward<_Mappedty>(_Mapval)));
		}

	template<class _Mappedty>
		iterator insert_or_assign(const_iterator, key_type&& _Keyval,
			_Mappedty&& _Mapval)
		{	
		return (_Insert_or_assign(::std:: move(_Keyval),
			::std:: forward<_Mappedty>(_Mapval)).first);
		}

	map(::std:: initializer_list<value_type> _Ilist)
		: _Mybase(key_compare())
		{	
		insert(_Ilist);
		}

	map(::std:: initializer_list<value_type> _Ilist,
		const key_compare& _Pred)
		: _Mybase(_Pred)
		{	
		insert(_Ilist);
		}

	map(::std:: initializer_list<value_type> _Ilist,
		const key_compare& _Pred, const allocator_type& _Al)
		: _Mybase(_Pred, _Al)
		{	
		insert(_Ilist);
		}

	_Myt& operator=(::std:: initializer_list<value_type> _Ilist)
		{	
		_Mybase::clear();
		insert(_Ilist);
		return (*this);
		}

	mapped_type& operator[](const key_type& _Keyval)
		{	
		return (try_emplace(_Keyval).first->second);
		}

	mapped_type& at(const key_type& _Keyval)
		{	
		iterator _Where = _Mybase::lower_bound(_Keyval);
		if (_Where == _Mybase::end()
			|| _Mybase::_Getcomp()(_Keyval, _Mybase::_Key(_Where._Mynode())))
			_Xout_of_range("invalid map<K, T> key");
		return (_Where->second);
		}

	const mapped_type& at(const key_type& _Keyval) const
		{	
		const_iterator _Where = _Mybase::lower_bound(_Keyval);
		if (_Where == _Mybase::end()
			|| _Mybase::_Getcomp()(_Keyval, _Mybase::_Key(_Where._Mynode())))
			_Xout_of_range("invalid map<K, T> key");
		return (_Where->second);
		}
	};

template<class _Kty,
	class _Ty,
	class _Pr,
	class _Alloc> inline
	void swap(map<_Kty, _Ty, _Pr, _Alloc>& _Left,
		map<_Kty, _Ty, _Pr, _Alloc>& _Right)
		noexcept(noexcept(_Left.swap(_Right)))
	{	
	_Left.swap(_Right);
	}

		
template<class _Kty,
	class _Ty,
	class _Pr = less<_Kty>,
	class _Alloc = allocator<pair<const _Kty, _Ty> > >
	class multimap
		: public _Tree<_Tmap_traits<_Kty, _Ty, _Pr, _Alloc, true> >
	{	
public:
	typedef multimap<_Kty, _Ty, _Pr, _Alloc> _Myt;
	typedef _Tree<_Tmap_traits<_Kty, _Ty, _Pr, _Alloc, true> > _Mybase;
	typedef _Kty key_type;
	typedef _Ty mapped_type;
	typedef _Pr key_compare;
	typedef typename _Mybase::value_compare value_compare;
	typedef typename _Mybase::allocator_type allocator_type;
	typedef typename _Mybase::size_type size_type;
	typedef typename _Mybase::difference_type difference_type;
	typedef typename _Mybase::pointer pointer;
	typedef typename _Mybase::const_pointer const_pointer;
	typedef typename _Mybase::reference reference;
	typedef typename _Mybase::const_reference const_reference;
	typedef typename _Mybase::iterator iterator;
	typedef typename _Mybase::const_iterator const_iterator;
	typedef typename _Mybase::reverse_iterator reverse_iterator;
	typedef typename _Mybase::const_reverse_iterator
		const_reverse_iterator;
	typedef typename _Mybase::value_type value_type;

	typedef typename _Mybase::_Alty _Alty;

	multimap()
		: _Mybase(key_compare())
		{	
		}

	explicit multimap(const allocator_type& _Al)
		: _Mybase(key_compare(), _Al)
		{	
		}

	multimap(const _Myt& _Right)
		: _Mybase(_Right,
			_Right._Getal().select_on_container_copy_construction())
		{	
		}

	multimap(const _Myt& _Right, const allocator_type& _Al)
		: _Mybase(_Right, _Al)
		{	
		}

	explicit multimap(const key_compare& _Pred)
		: _Mybase(_Pred)
		{	
		}

	multimap(const key_compare& _Pred, const allocator_type& _Al)
		: _Mybase(_Pred, _Al)
		{	
		}

	template<class _Iter>
		multimap(_Iter _First, _Iter _Last)
		: _Mybase(key_compare())
		{	
		insert(_First, _Last);
		}

	template<class _Iter>
		multimap(_Iter _First, _Iter _Last,
			const key_compare& _Pred)
		: _Mybase(_Pred)
		{	
		insert(_First, _Last);
		}

	template<class _Iter>
		multimap(_Iter _First, _Iter _Last,
			const key_compare& _Pred, const allocator_type& _Al)
		: _Mybase(_Pred, _Al)
		{	
		insert(_First, _Last);
		}

	_Myt& operator=(const _Myt& _Right)
		{	
		_Mybase::operator=(_Right);
		return (*this);
		}

	multimap(_Myt&& _Right)
		: _Mybase(::std:: move(_Right))
		{	
		}

	multimap(_Myt&& _Right, const allocator_type& _Al)
		: _Mybase(::std:: move(_Right), _Al)
		{	
		}

	_Myt& operator=(_Myt&& _Right)
		noexcept(_Alty::is_always_equal::value && is_nothrow_move_assignable<_Pr>::value)
		{	
		_Mybase::operator=(::std:: move(_Right));
		return (*this);
		}

	template<class... _Valty>
		iterator emplace(_Valty&&... _Val)
		{	
		return (_Mybase::emplace(::std:: forward<_Valty>(_Val)...).first);
		}

	void swap(_Myt& _Right)
		noexcept(_Alty::is_always_equal::value && _Is_nothrow_swappable<_Pr>::value)
		{	
		_Mybase::swap(_Right);
		}

	using _Mybase::insert;

	template<class _Valty,
		class = enable_if_t<is_constructible<value_type, _Valty>::value> >
		iterator insert(_Valty&& _Val)
		{	
		return (this->emplace(::std:: forward<_Valty>(_Val)));
		}

	template<class _Valty,
		class = enable_if_t<is_constructible<value_type, _Valty>::value> >
		iterator insert(const_iterator _Where, _Valty&& _Val)
		{	
		return (this->emplace_hint(_Where, ::std:: forward<_Valty>(_Val)));
		}

	multimap(::std:: initializer_list<value_type> _Ilist)
		: _Mybase(key_compare())
		{	
		insert(_Ilist);
		}

	multimap(::std:: initializer_list<value_type> _Ilist,
		const key_compare& _Pred)
		: _Mybase(_Pred)
		{	
		insert(_Ilist);
		}

	multimap(::std:: initializer_list<value_type> _Ilist,
		const key_compare& _Pred, const allocator_type& _Al)
		: _Mybase(_Pred, _Al)
		{	
		insert(_Ilist);
		}

	_Myt& operator=(::std:: initializer_list<value_type> _Ilist)
		{	
		_Mybase::clear();
		insert(_Ilist);
		return (*this);
		}
	};

template<class _Kty,
	class _Ty,
	class _Pr,
	class _Alloc> inline
	void swap(multimap<_Kty, _Ty, _Pr, _Alloc>& _Left,
		multimap<_Kty, _Ty, _Pr, _Alloc>& _Right)
		noexcept(noexcept(_Left.swap(_Right)))
	{	
	_Left.swap(_Right);
	}
}
 
 #pragma warning(pop)
 #pragma pack(pop)










#pragma once





 #pragma pack(push,8)
 #pragma warning(push,3)
 
 
 #pragma warning(disable: 4244 28309 28285)

namespace std {
		
const int _ISORT_MAX = 32;	

template<class _Iter1,
	class _Iter2,
	class _UIter1,
	class _UIter2> inline
	pair<_Iter1, _Iter2>
		_Rechecked_both(_Iter1 _Dest1, _Iter2 _Dest2, pair<_UIter1, _UIter2> _Src)
	{	
	return (pair<_Iter1, _Iter2>(
		_Rechecked(_Dest1, _Src.first),
		_Rechecked(_Dest2, _Src.second)
		));
	}

 































		
template<class _InIt,
	class _Fn1> inline
	void _For_each_unchecked(_InIt _First, _InIt _Last, _Fn1& _Func)
	{	
	for (; _First != _Last; ++_First)
		_Func(*_First);
	}

template<class _InIt,
	class _Fn1> inline
	_Fn1 for_each(_InIt _First, _InIt _Last, _Fn1 _Func)
	{	
	;
	_For_each_unchecked(_Unchecked(_First), _Unchecked(_Last), _Func);
	return (_Func);
	}

		
template<class _InIt,
	class _Pr> inline
	_InIt _Find_if_unchecked(_InIt _First, _InIt _Last, _Pr& _Pred)
	{	
	for (; _First != _Last; ++_First)
		if (_Pred(*_First))
			break;
	return (_First);
	}

template<class _InIt,
	class _Pr> inline
	_InIt find_if(_InIt _First, _InIt _Last, _Pr _Pred)
	{	
	;
	return (_Rechecked(_First,
		_Find_if_unchecked(_Unchecked(_First), _Unchecked(_Last), _Pred)));
	}

		
template<class _FwdIt,
	class _Pr> inline
	_FwdIt _Adjacent_find_unchecked(_FwdIt _First, _FwdIt _Last, _Pr& _Pred)
	{	
	if (_First != _Last)
		for (_FwdIt _Firstb; (void)(_Firstb = _First), ++_First != _Last; )
			if (_Pred(*_Firstb, *_First))
				return (_Firstb);
	return (_Last);
	}

template<class _FwdIt,
	class _Pr> inline
	_FwdIt adjacent_find(_FwdIt _First, _FwdIt _Last, _Pr _Pred)
	{	
	;
	;
	return (_Rechecked(_First,
		_Adjacent_find_unchecked(_Unchecked(_First), _Unchecked(_Last), _Pred)));
	}

		
template<class _FwdIt> inline
	_FwdIt adjacent_find(_FwdIt _First, _FwdIt _Last)
	{	
	return (::std:: adjacent_find(_First, _Last, equal_to<>()));
	}

		
template<class _InIt,
	class _Pr> inline
	typename iterator_traits<_InIt>::difference_type
		_Count_if_unchecked(_InIt _First, _InIt _Last, _Pr& _Pred)
	{	
	typename iterator_traits<_InIt>::difference_type _Count = 0;

	for (; _First != _Last; ++_First)
		if (_Pred(*_First))
			++_Count;
	return (_Count);
	}

template<class _InIt,
	class _Pr> inline
	typename iterator_traits<_InIt>::difference_type
		count_if(_InIt _First, _InIt _Last, _Pr _Pred)
	{	
	;
	return (_Count_if_unchecked(_Unchecked(_First), _Unchecked(_Last), _Pred));
	}

		
template<class _InIt1,
	class _InIt2,
	class _Pr> inline
	pair<_InIt1, _InIt2>
		_Mismatch_unchecked(_InIt1 _First1, _InIt1 _Last1,
			_InIt2 _First2, _Pr& _Pred)
	{	
	for (; _First1 != _Last1 && _Pred(*_First1, *_First2); )
		{	
		++_First1;
		++_First2;
		}

	return (pair<_InIt1, _InIt2>(_First1, _First2));
	}

template<class _InIt1,
	class _InIt2,
	class _Pr> inline
	auto _Mismatch_no_deprecate1(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _Pr& _Pred, input_iterator_tag, input_iterator_tag)
			-> pair<_InIt1, decltype(_Unchecked_idl0(_First2))>
	{	
	return (_Mismatch_unchecked(_First1, _Last1,
		_Unchecked_idl0(_First2), _Pred));
	}

template<class _InIt1,
	class _InIt2,
	class _Pr> inline
	pair<_InIt1, decltype(_Unchecked(::std:: declval<_InIt2>()))>
		_Mismatch_no_deprecate1(_InIt1 _First1, _InIt1 _Last1,
			_InIt2 _First2, _Pr& _Pred, random_access_iterator_tag, random_access_iterator_tag)
	{	
	;
	return (_Mismatch_unchecked(_First1, _Last1, _Unchecked(_First2), _Pred));
	}

template<class _InIt1,
	class _InIt2,
	class _Pr> inline
	pair<_InIt1, _InIt2>
		_Mismatch_no_deprecate(_InIt1 _First1, _InIt1 _Last1,
			_InIt2 _First2, _Pr& _Pred)
	{	
	;
	;
	return (_Rechecked_both(_First1, _First2,
		_Mismatch_no_deprecate1(_Unchecked(_First1), _Unchecked(_Last1),
			_First2, _Pred, _Iter_cat_t<_InIt1>(), _Iter_cat_t<_InIt2>())));
	}

template<class _InIt1,
	class _InIt2,
	class _Pr> inline
	pair<_InIt1, _InIt2>
		mismatch(_InIt1 _First1, _InIt1 _Last1,
			_InIt2 _First2, _Pr _Pred)
	{	
	struct _Unchecked_iterators { static void  _Deprecate(false_type) { } static void _Deprecate(true_type) { } }; (_Unchecked_iterators::_Deprecate(_Is_checked(_First2)));
	return (_Mismatch_no_deprecate(_First1, _Last1, _First2, _Pred));
	}

 














		
template<class _InIt1,
	class _InIt2> inline
	pair<_InIt1, _InIt2>
		mismatch(_InIt1 _First1, _InIt1 _Last1,
			_InIt2 _First2)
	{	
	return (::std:: mismatch(_First1, _Last1, _First2,
		equal_to<>()));
	}

 












		
template<class _InIt1,
	class _InIt2,
	class _Pr> inline
	pair<_InIt1, _InIt2>
		_Mismatch_unchecked(_InIt1 _First1, _InIt1 _Last1,
			_InIt2 _First2, _InIt2 _Last2, _Pr& _Pred)
	{	
	for (; _First1 != _Last1 && _First2 != _Last2
		&& _Pred(*_First1, *_First2); )
		{	
		++_First1;
		++_First2;
		}

	return (pair<_InIt1, _InIt2>(_First1, _First2));
	}

template<class _InIt1,
	class _InIt2,
	class _Pr> inline
	pair<_InIt1, _InIt2>
		mismatch(_InIt1 _First1, _InIt1 _Last1,
			_InIt2 _First2, _InIt2 _Last2, _Pr _Pred)
	{	
	;
	;
	;
	return (_Rechecked_both(_First1, _First2,
		_Mismatch_unchecked(_Unchecked(_First1), _Unchecked(_Last1),
			_Unchecked(_First2), _Unchecked(_Last2), _Pred)));
	}

		
template<class _InIt1,
	class _InIt2> inline
	pair<_InIt1, _InIt2>
		mismatch(_InIt1 _First1, _InIt1 _Last1,
			_InIt2 _First2, _InIt2 _Last2)
	{	
	return (::std:: mismatch(_First1, _Last1, _First2, _Last2,
		equal_to<>()));
	}

		
template<class _InIt,
	class _Pr> inline
	bool _All_of_unchecked(_InIt _First, _InIt _Last, _Pr& _Pred)
	{	
	for (; _First != _Last; ++_First)
		if (!_Pred(*_First))
			return (false);
	return (true);
	}

template<class _InIt,
	class _Pr> inline
	bool all_of(_InIt _First, _InIt _Last, _Pr _Pred)
	{	
	;
	return (_All_of_unchecked(_Unchecked(_First), _Unchecked(_Last), _Pred));
	}

		
template<class _InIt,
	class _Pr> inline
	bool _Any_of_unchecked(_InIt _First, _InIt _Last, _Pr& _Pred)
	{	
	for (; _First != _Last; ++_First)
		if (_Pred(*_First))
			return (true);
	return (false);
	}

template<class _InIt,
	class _Pr> inline
	bool any_of(_InIt _First, _InIt _Last, _Pr _Pred)
	{	
	;
	return (_Any_of_unchecked(_Unchecked(_First), _Unchecked(_Last), _Pred));
	}

		
template<class _InIt,
	class _Pr> inline
	bool _None_of_unchecked(_InIt _First, _InIt _Last, _Pr& _Pred)
	{	
	for (; _First != _Last; ++_First)
		if (_Pred(*_First))
			return (false);
	return (true);
	}

template<class _InIt,
	class _Pr> inline
	bool none_of(_InIt _First, _InIt _Last, _Pr _Pred)
	{	
	;
	return (_None_of_unchecked(_Unchecked(_First), _Unchecked(_Last), _Pred));
	}

		
template<class _InIt,
	class _Pr> inline
	_InIt _Find_if_not_unchecked(_InIt _First, _InIt _Last, _Pr& _Pred)
	{	
	for (; _First != _Last; ++_First)
		if (!_Pred(*_First))
			break;
	return (_First);
	}

template<class _InIt,
	class _Pr> inline
	_InIt find_if_not(_InIt _First, _InIt _Last, _Pr _Pred)
	{	
	;
	return (_Rechecked(_First,
		_Find_if_not_unchecked(_Unchecked(_First), _Unchecked(_Last), _Pred)));
	}

		
template<class _InIt,
	class _OutIt,
	class _Pr> inline
	_OutIt _Copy_if_unchecked(_InIt _First, _InIt _Last, _OutIt _Dest,
		_Pr& _Pred)
	{	
	for (; _First != _Last; ++_First)
		if (_Pred(*_First))
			{	
			;
			*_Dest++ = *_First;
			}

	return (_Dest);
	}

template<class _InIt,
	class _OutIt,
	class _Pr> inline
	_OutIt _Copy_if_no_deprecate(_InIt _First, _InIt _Last, _OutIt _Dest,
		_Pr& _Pred)
	{	
	;
	return (_Rechecked(_Dest,
		_Copy_if_unchecked(_Unchecked(_First), _Unchecked(_Last),
		_Unchecked_idl0(_Dest), _Pred)));
	}

template<class _InIt,
	class _OutIt,
	class _Pr> inline
	_OutIt copy_if(_InIt _First, _InIt _Last, _OutIt _Dest,
		_Pr _Pred)
	{	
	struct _Unchecked_iterators { static void  _Deprecate(false_type) { } static void _Deprecate(true_type) { } }; (_Unchecked_iterators::_Deprecate(_Is_checked(_Dest)));
	return (_Copy_if_no_deprecate(_First, _Last, _Dest, _Pred));
	}

 













		
template<class _InIt,
	class _OutIt1,
	class _OutIt2,
	class _Pr> inline
	pair<_OutIt1, _OutIt2>
		_Partition_copy_unchecked(_InIt _First, _InIt _Last,
			_OutIt1 _Dest1, _OutIt2 _Dest2, _Pr& _Pred)
	{	
	for (; _First != _Last; ++_First)
		if (_Pred(*_First))
			{	
			;
			*_Dest1++ = *_First;
			}
		else
			{	
			;
			*_Dest2++ = *_First;
			}

	return (pair<_OutIt1, _OutIt2>(_Dest1, _Dest2));
	}

template<class _InIt,
	class _OutIt1,
	class _OutIt2,
	class _Pr> inline
	pair<_OutIt1, _OutIt2>
		_Partition_copy_no_deprecate(_InIt _First, _InIt _Last,
			_OutIt1 _Dest1, _OutIt2 _Dest2, _Pr& _Pred)
	{	
	;
	return (_Rechecked_both(_Dest1, _Dest2,
		_Partition_copy_unchecked(_Unchecked(_First), _Unchecked(_Last),
		_Unchecked_idl0(_Dest1), _Unchecked_idl0(_Dest2), _Pred)));
	}

template<class _InIt,
	class _OutIt1,
	class _OutIt2,
	class _Pr> inline
	pair<_OutIt1, _OutIt2>
		partition_copy(_InIt _First, _InIt _Last,
			_OutIt1 _Dest1, _OutIt2 _Dest2, _Pr _Pred)
	{	
	struct _Unchecked_iterators { static void  _Deprecate(false_type) { } static void _Deprecate(true_type) { } };
	(_Unchecked_iterators::_Deprecate(_Is_checked(_Dest1)));
	(_Unchecked_iterators::_Deprecate(_Is_checked(_Dest2)));
	return (_Partition_copy_no_deprecate(_First, _Last, _Dest1, _Dest2, _Pred));
	}

 
















































		
template<class _InIt,
	class _Pr> inline
	bool _Is_partitioned_unchecked(_InIt _First, _InIt _Last, _Pr& _Pred)
	{	
	for (; _First != _Last; ++_First)
		if (!_Pred(*_First))
			break;	
	for (; _First != _Last; ++_First)
		if (_Pred(*_First))
			return (false);	
	return (true);
	}

template<class _InIt,
	class _Pr> inline
	bool is_partitioned(_InIt _First, _InIt _Last, _Pr _Pred)
	{	
	;
	return (_Is_partitioned_unchecked(_Unchecked(_First), _Unchecked(_Last),
		_Pred));
	}

		
template<class _FwdIt,
	class _Pr> inline
	_FwdIt _Partition_point_unchecked(_FwdIt _First, _FwdIt _Last, _Pr& _Pred)
	{	
	_Iter_diff_t<_FwdIt> _Count = ::std:: distance(_First, _Last);
	while (0 < _Count)
		{	
		_Iter_diff_t<_FwdIt> _Count2 = _Count / 2;
		_FwdIt _Mid = _First;
		::std:: advance(_Mid, _Count2);

		if (_Pred(*_Mid))
			{	
			_First = ++_Mid;
			_Count -= _Count2 + 1;
			}
		else
			_Count = _Count2;
		}

	return (_First);
	}

template<class _FwdIt,
	class _Pr> inline
	_FwdIt partition_point(_FwdIt _First, _FwdIt _Last, _Pr _Pred)
	{	
	;
	return (_Rechecked(_First,
		_Partition_point_unchecked(_Unchecked(_First), _Unchecked(_Last), _Pred)));
	}

		
template<class _FwdIt1,
	class _FwdIt2,
	class _Pr> inline
	_FwdIt1 _Search_unchecked(_FwdIt1 _First1, _FwdIt1 _Last1,
		_FwdIt2 _First2, _FwdIt2 _Last2, _Pr& _Pred,
		forward_iterator_tag, forward_iterator_tag)
	{	
	for (; ; ++_First1)
		{	
		_FwdIt1 _Mid1 = _First1;
		for (_FwdIt2 _Mid2 = _First2; ; ++_Mid1, (void)++_Mid2)
			if (_Mid2 == _Last2)
				return (_First1);
			else if (_Mid1 == _Last1)
				return (_Last1);
			else if (!_Pred(*_Mid1, *_Mid2))
				break;
		}
	}

template<class _FwdIt1,
	class _FwdIt2,
	class _Pr> inline
	_FwdIt1 _Search_unchecked(_FwdIt1 _First1, _FwdIt1 _Last1,
		_FwdIt2 _First2, _FwdIt2 _Last2, _Pr& _Pred,
		random_access_iterator_tag, random_access_iterator_tag)
	{	
	_Iter_diff_t<_FwdIt1> _Count1 = _Last1 - _First1;
	_Iter_diff_t<_FwdIt2> _Count2 = _Last2 - _First2;

	for (; _Count2 <= _Count1; ++_First1, (void)--_Count1)
		{	
		_FwdIt1 _Mid1 = _First1;
		for (_FwdIt2 _Mid2 = _First2; ; ++_Mid1, (void)++_Mid2)
			if (_Mid2 == _Last2)
				return (_First1);
			else if (!_Pred(*_Mid1, *_Mid2))
				break;
		}

	return (_Last1);
	}

template<class _FwdIt1,
	class _FwdIt2,
	class _Pr> inline
	_FwdIt1 search(_FwdIt1 _First1, _FwdIt1 _Last1,
		_FwdIt2 _First2, _FwdIt2 _Last2, _Pr _Pred)
	{	
	;
	;
	;
	return (_Rechecked(_First1,
		_Search_unchecked(_Unchecked(_First1), _Unchecked(_Last1),
			_Unchecked(_First2), _Unchecked(_Last2), _Pred,
			_Iter_cat_t<_FwdIt1>(), _Iter_cat_t<_FwdIt2>())));
	}

		
template<class _FwdIt1,
	class _FwdIt2> inline
	_FwdIt1 search(_FwdIt1 _First1, _FwdIt1 _Last1,
		_FwdIt2 _First2, _FwdIt2 _Last2)
	{	
	return (::std:: search(_First1, _Last1, _First2, _Last2,
		equal_to<>()));
	}

		
template<class _FwdIt,
	class _Diff,
	class _Ty,
	class _Pr> inline
	_FwdIt _Search_n_unchecked(_FwdIt _First, _FwdIt _Last,
		_Diff _Count, const _Ty& _Val, _Pr& _Pred, forward_iterator_tag)
	{	
	if (_Count <= 0)
		return (_First);

	for (; _First != _Last; ++_First)
		if (_Pred(*_First, _Val))
			{	
			_FwdIt _Mid = _First;

			for (_Diff _Count1 = _Count; ; )
				if (--_Count1 == 0)
					return (_First);	
				else if (++_Mid == _Last)
					return (_Last);	
				else if (!_Pred(*_Mid, _Val))
					{	
					break;
					}

			_First = _Mid;	
			}

	return (_Last);
	}

template<class _FwdIt,
	class _Diff,
	class _Ty,
	class _Pr> inline
	_FwdIt _Search_n_unchecked(_FwdIt _First, _FwdIt _Last,
		_Diff _Count, const _Ty& _Val, _Pr& _Pred, random_access_iterator_tag)
	{	
	if (_Count <= 0)
		return (_First);

	_FwdIt _Oldfirst = _First;
	for (_Diff _Inc = 0; _Count <= _Last - _Oldfirst; )
		{	
		_First = _Oldfirst + _Inc;
		if (_Pred(*_First, _Val))
			{	
			_Diff _Count1 = _Count;
			_FwdIt _Mid = _First;

			for (; _Oldfirst != _First && _Pred(_First[-1], _Val);
				--_First)
				--_Count1;	

			if (_Count1 <= _Last - _Mid)
				for (; ; )
					{	
					if (--_Count1 == 0)
						return (_First);	
					else if (!_Pred(*++_Mid, _Val))
						{	
						break;
						}
					}
			_Oldfirst = ++_Mid;	
			_Inc = 0;
			}
		else
			{	
			_Oldfirst = _First + 1;
			_Inc = _Count - 1;
			}
		}

	return (_Last);
	}

template<class _FwdIt,
	class _Diff,
	class _Ty,
	class _Pr> inline
	_FwdIt search_n(_FwdIt _First, _FwdIt _Last,
		_Diff _Count, const _Ty& _Val, _Pr _Pred)
	{	
	;
	return (_Rechecked(_First,
		_Search_n_unchecked(_Unchecked(_First), _Unchecked(_Last), _Count, _Val,
			_Pred, _Iter_cat_t<_FwdIt>())));
	}

		
template<class _FwdIt,
	class _Diff,
	class _Ty> inline
	_FwdIt search_n(_FwdIt _First, _FwdIt _Last,
		_Diff _Count, const _Ty& _Val)
	{	
	return (::std:: search_n(_First, _Last, _Count, _Val,
		equal_to<>()));
	}

		
template<class _FwdIt1,
	class _FwdIt2,
	class _Pr> inline
	_FwdIt1 _Find_end_unchecked(_FwdIt1 _First1, _FwdIt1 _Last1,
		_FwdIt2 _First2, _FwdIt2 _Last2, _Pr& _Pred)
	{	
	_Iter_diff_t<_FwdIt1> _Count1 = ::std:: distance(_First1, _Last1);
	_Iter_diff_t<_FwdIt2> _Count2 = ::std:: distance(_First2, _Last2);
	_FwdIt1 _Ans = _Last1;

	if (0 < _Count2)
		{	
		;
		for (; _Count2 <= _Count1; ++_First1, (void)--_Count1)
			{	
			_FwdIt1 _Mid1 = _First1;
			for (_FwdIt2 _Mid2 = _First2; ; ++_Mid1)
				if (!_Pred(*_Mid1, *_Mid2))
					break;
				else if (++_Mid2 == _Last2)
					{	
					_Ans = _First1;
					break;
					}
			}
		}

	return (_Ans);
	}

template<class _FwdIt1,
	class _FwdIt2,
	class _Pr> inline
	_FwdIt1 find_end(_FwdIt1 _First1, _FwdIt1 _Last1,
		_FwdIt2 _First2, _FwdIt2 _Last2, _Pr _Pred)
	{	
	;
	;
	return (_Rechecked(_First1,
		_Find_end_unchecked(_Unchecked(_First1), _Unchecked(_Last1),
			_Unchecked(_First2), _Unchecked(_Last2), _Pred)));
	}

		
template<class _FwdIt1,
	class _FwdIt2> inline
	_FwdIt1 find_end(_FwdIt1 _First1, _FwdIt1 _Last1,
		_FwdIt2 _First2, _FwdIt2 _Last2)
	{	
	return (::std:: find_end(_First1, _Last1, _First2, _Last2,
		equal_to<>()));
	}

		
template<class _FwdIt1,
	class _FwdIt2,
	class _Pr> inline
	_FwdIt1 _Find_first_of_unchecked(_FwdIt1 _First1, _FwdIt1 _Last1,
		_FwdIt2 _First2, _FwdIt2 _Last2, _Pr& _Pred)
	{	
	for (; _First1 != _Last1; ++_First1)
		for (_FwdIt2 _Mid2 = _First2; _Mid2 != _Last2; ++_Mid2)
			if (_Pred(*_First1, *_Mid2))
				return (_First1);
	return (_First1);
	}

template<class _FwdIt1,
	class _FwdIt2,
	class _Pr> inline
	_FwdIt1 find_first_of(_FwdIt1 _First1, _FwdIt1 _Last1,
		_FwdIt2 _First2, _FwdIt2 _Last2, _Pr _Pred)
	{	
	;
	;
	;
	return (_Rechecked(_First1,
		_Find_first_of_unchecked(_Unchecked(_First1), _Unchecked(_Last1),
			_Unchecked(_First2), _Unchecked(_Last2), _Pred)));
	}

		
template<class _FwdIt1,
	class _FwdIt2> inline
	_FwdIt1 find_first_of(_FwdIt1 _First1, _FwdIt1 _Last1,
		_FwdIt2 _First2, _FwdIt2 _Last2)
	{	
	return (::std:: find_first_of(_First1, _Last1, _First2, _Last2,
		equal_to<>()));
	}

		
template<class _FwdIt1,
	class _FwdIt2> inline
	_FwdIt2 _Swap_ranges_unchecked(_FwdIt1 _First1, _FwdIt1 _Last1,
		_FwdIt2 _Dest)
	{	
	for (; _First1 != _Last1; ++_First1, (void)++_Dest)
		::std:: iter_swap(_First1, _Dest);
	return (_Dest);
	}

template<class _FwdIt1,
	class _FwdIt2> inline
	_FwdIt2 _Swap_ranges1(_FwdIt1 _First1, _FwdIt1 _Last1,
		_FwdIt2 _Dest,
		forward_iterator_tag, forward_iterator_tag)
	{	
	return (_Rechecked(_Dest,
		_Swap_ranges_unchecked(_First1, _Last1, _Unchecked_idl0(_Dest))));
	}

template<class _FwdIt1,
	class _FwdIt2> inline
	_FwdIt2 _Swap_ranges1(_FwdIt1 _First1, _FwdIt1 _Last1,
		_FwdIt2 _Dest,
		random_access_iterator_tag, random_access_iterator_tag)
	{	
	;
	return (_Rechecked(_Dest,
		_Swap_ranges_unchecked(_First1, _Last1, _Unchecked(_Dest))));
	}

template<class _FwdIt1,
	class _FwdIt2> inline
	_FwdIt2 swap_ranges(_FwdIt1 _First1, _FwdIt1 _Last1,
		_FwdIt2 _Dest)
	{	
	struct _Unchecked_iterators { static void  _Deprecate(false_type) { } static void _Deprecate(true_type) { } }; (_Unchecked_iterators::_Deprecate(_Is_checked(_Dest)));
	;
	return (_Swap_ranges1(_Unchecked(_First1), _Unchecked(_Last1),
		_Dest, _Iter_cat_t<_FwdIt1>(), _Iter_cat_t<_FwdIt2>()));
	}

 












		
template<class _InIt,
	class _OutIt,
	class _Fn1> inline
	_OutIt _Transform_unchecked(_InIt _First, _InIt _Last,
		_OutIt _Dest, _Fn1& _Func)
	{	
	for (; _First != _Last; ++_First, (void)++_Dest)
		*_Dest = _Func(*_First);
	return (_Dest);
	}

template<class _InIt,
	class _OutIt,
	class _Fn1> inline
	_OutIt _Transform_no_deprecate1(_InIt _First, _InIt _Last,
		_OutIt _Dest, _Fn1& _Func,
		input_iterator_tag, _Any_tag)
	{	
	return (_Rechecked(_Dest,
		_Transform_unchecked(_First, _Last, _Unchecked_idl0(_Dest), _Func)));
	}

template<class _InIt,
	class _OutIt,
	class _Fn1> inline
	_OutIt _Transform_no_deprecate1(_InIt _First, _InIt _Last,
		_OutIt _Dest, _Fn1& _Func,
		random_access_iterator_tag, random_access_iterator_tag)
	{	
	;
	return (_Rechecked(_Dest,
		_Transform_unchecked(_First, _Last, _Unchecked(_Dest), _Func)));
	}

template<class _InIt,
	class _OutIt,
	class _Fn1> inline
	_OutIt _Transform_no_deprecate(_InIt _First, _InIt _Last,
		_OutIt _Dest, _Fn1& _Func)
	{	
	;
	;
	return (_Transform_no_deprecate1(_Unchecked(_First), _Unchecked(_Last),
		_Dest, _Func, _Iter_cat_t<_InIt>(), _Iter_cat_t<_OutIt>()));
	}

template<class _InIt,
	class _OutIt,
	class _Fn1> inline
	_OutIt transform(_InIt _First, _InIt _Last,
		_OutIt _Dest, _Fn1 _Func)
	{	
	struct _Unchecked_iterators { static void  _Deprecate(false_type) { } static void _Deprecate(true_type) { } }; (_Unchecked_iterators::_Deprecate(_Is_checked(_Dest)));
	return (_Transform_no_deprecate(_First, _Last, _Dest, _Func));
	}

 













		
template<class _InIt1,
	class _InIt2,
	class _OutIt,
	class _Fn2> inline
	_OutIt _Transform_unchecked(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _OutIt _Dest, _Fn2& _Func)
	{	
	for (; _First1 != _Last1; ++_First1, (void)++_First2, ++_Dest)
		*_Dest = _Func(*_First1, *_First2);
	return (_Dest);
	}

template<class _InIt1,
	class _InIt2,
	class _OutIt,
	class _Fn2> inline
	_OutIt _Transform_no_deprecate2(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _OutIt _Dest, _Fn2& _Func,
		input_iterator_tag, _Any_tag)
	{	
	return (_Rechecked(_Dest,
		_Transform_unchecked(_First1, _Last1, _First2, _Unchecked_idl0(_Dest), _Func)));
	}

template<class _InIt1,
	class _InIt2,
	class _OutIt,
	class _Fn2> inline
	_OutIt _Transform_no_deprecate2(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _OutIt _Dest, _Fn2& _Func,
		random_access_iterator_tag, random_access_iterator_tag)
	{	
	;
	return (_Rechecked(_Dest,
		_Transform_unchecked(_First1, _Last1, _First2, _Unchecked(_Dest), _Func)));
	}

template<class _InIt1,
	class _InIt2,
	class _OutIt,
	class _Fn2> inline
	_OutIt _Transform_no_deprecate1(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _OutIt _Dest, _Fn2& _Func,
		input_iterator_tag, input_iterator_tag)
	{	
	return (_Transform_no_deprecate2(_First1, _Last1,
		_Unchecked_idl0(_First2), _Dest, _Func,
		_Iter_cat_t<_InIt1>(), _Iter_cat_t<_OutIt>()));
	}

template<class _InIt1,
	class _InIt2,
	class _OutIt,
	class _Fn2> inline
	_OutIt _Transform_no_deprecate1(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _OutIt _Dest, _Fn2& _Func,
		random_access_iterator_tag, random_access_iterator_tag)
	{	
	;
	return (_Transform_no_deprecate2(_First1, _Last1,
		_Unchecked(_First2), _Dest, _Func,
		_Iter_cat_t<_InIt1>(), _Iter_cat_t<_OutIt>()));
	}

template<class _InIt1,
	class _InIt2,
	class _OutIt,
	class _Fn2> inline
	_OutIt _Transform_no_deprecate(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _OutIt _Dest, _Fn2& _Func)
	{	
	;
	;
	;
	return (_Transform_no_deprecate1(_Unchecked(_First1), _Unchecked(_Last1),
		_First2, _Dest, _Func,
		_Iter_cat_t<_InIt1>(), _Iter_cat_t<_InIt2>()));
	}

template<class _InIt1,
	class _InIt2,
	class _OutIt,
	class _Fn2> inline
	_OutIt transform(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _OutIt _Dest, _Fn2 _Func)
	{	
	struct _Unchecked_iterators { static void  _Deprecate(false_type) { } static void _Deprecate(true_type) { } };
	(_Unchecked_iterators::_Deprecate(_Is_checked(_First2)));
	(_Unchecked_iterators::_Deprecate(_Is_checked(_Dest)));
	return (_Transform_no_deprecate(_First1, _Last1, _First2, _Dest, _Func));
	}

 











































		
template<class _FwdIt,
	class _Ty> inline
	void _Replace_unchecked(_FwdIt _First, _FwdIt _Last,
		const _Ty& _Oldval, const _Ty& _Newval)
	{	
	for (; _First != _Last; ++_First)
		if (*_First == _Oldval)
			*_First = _Newval;
	}

template<class _FwdIt,
	class _Ty> inline
	void replace(_FwdIt _First, _FwdIt _Last,
		const _Ty& _Oldval, const _Ty& _Newval)
	{	
	;
	_Replace_unchecked(_Unchecked(_First), _Unchecked(_Last),
		_Oldval, _Newval);
	}

		
template<class _FwdIt,
	class _Pr,
	class _Ty> inline
	void _Replace_if_unchecked(_FwdIt _First, _FwdIt _Last, _Pr& _Pred, const _Ty& _Val)
	{	
	for (; _First != _Last; ++_First)
		if (_Pred(*_First))
			*_First = _Val;
	}

template<class _FwdIt,
	class _Pr,
	class _Ty> inline
	void replace_if(_FwdIt _First, _FwdIt _Last, _Pr _Pred, const _Ty& _Val)
	{	
	;
	_Replace_if_unchecked(_Unchecked(_First), _Unchecked(_Last),
		_Pred, _Val);
	}

		
template<class _InIt,
	class _OutIt,
	class _Ty> inline
	_OutIt _Replace_copy_unchecked(_InIt _First, _InIt _Last,
		_OutIt _Dest, const _Ty& _Oldval, const _Ty& _Newval)
	{	
	for (; _First != _Last; ++_First, (void)++_Dest)
		*_Dest = *_First == _Oldval ? _Newval : *_First;
	return (_Dest);
	}

template<class _InIt,
	class _OutIt,
	class _Ty> inline
	_OutIt _Replace_copy1(_InIt _First, _InIt _Last,
		_OutIt _Dest, const _Ty& _Oldval, const _Ty& _Newval,
		input_iterator_tag, _Any_tag)
	{	
	return (_Rechecked(_Dest,
		_Replace_copy_unchecked(_First, _Last, _Unchecked_idl0(_Dest),
		_Oldval, _Newval)));
	}

template<class _InIt,
	class _OutIt,
	class _Ty> inline
	_OutIt _Replace_copy1(_InIt _First, _InIt _Last,
		_OutIt _Dest, const _Ty& _Oldval, const _Ty& _Newval,
		random_access_iterator_tag, random_access_iterator_tag)
	{	
	;
	return (_Rechecked(_Dest,
		_Replace_copy_unchecked(_First, _Last, _Unchecked(_Dest),
		_Oldval, _Newval)));
	}

template<class _InIt,
	class _OutIt,
	class _Ty> inline
	_OutIt replace_copy(_InIt _First, _InIt _Last,
		_OutIt _Dest, const _Ty& _Oldval, const _Ty& _Newval)
	{	
	struct _Unchecked_iterators { static void  _Deprecate(false_type) { } static void _Deprecate(true_type) { } }; (_Unchecked_iterators::_Deprecate(_Is_checked(_Dest)));
	;
	return (_Replace_copy1(_Unchecked(_First), _Unchecked(_Last),
		_Dest, _Oldval, _Newval,
		_Iter_cat_t<_InIt>(), _Iter_cat_t<_OutIt>()));
	}

 













		
template<class _InIt,
	class _OutIt,
	class _Pr,
	class _Ty> inline
	_OutIt _Replace_copy_if_unchecked(_InIt _First, _InIt _Last,
		_OutIt _Dest, _Pr& _Pred, const _Ty& _Val)
	{	
	for (; _First != _Last; ++_First, (void)++_Dest)
		*_Dest = _Pred(*_First) ? _Val : *_First;
	return (_Dest);
	}

template<class _InIt,
	class _OutIt,
	class _Pr,
	class _Ty> inline
	_OutIt _Replace_copy_if_no_deprecate1(_InIt _First, _InIt _Last,
		_OutIt _Dest, _Pr& _Pred, const _Ty& _Val,
		input_iterator_tag, _Any_tag)
	{	
	return (_Rechecked(_Dest,
		_Replace_copy_if_unchecked(_First, _Last, _Unchecked_idl0(_Dest), _Pred, _Val)));
	}

template<class _InIt,
	class _OutIt,
	class _Pr,
	class _Ty> inline
	_OutIt _Replace_copy_if_no_deprecate1(_InIt _First, _InIt _Last,
		_OutIt _Dest, _Pr& _Pred, const _Ty& _Val,
		random_access_iterator_tag, random_access_iterator_tag)
	{	
	;
	return (_Rechecked(_Dest,
		_Replace_copy_if_unchecked(_First, _Last, _Unchecked(_Dest), _Pred, _Val)));
	}

template<class _InIt,
	class _OutIt,
	class _Pr,
	class _Ty> inline
	_OutIt _Replace_copy_if_no_deprecate(_InIt _First, _InIt _Last,
		_OutIt _Dest, _Pr& _Pred, const _Ty& _Val)
	{	
	;
	;
	return (_Replace_copy_if_no_deprecate1(_Unchecked(_First), _Unchecked(_Last),
		_Dest, _Pred, _Val,
		_Iter_cat_t<_InIt>(), _Iter_cat_t<_OutIt>()));
	}

template<class _InIt,
	class _OutIt,
	class _Pr,
	class _Ty> inline
	_OutIt replace_copy_if(_InIt _First, _InIt _Last,
		_OutIt _Dest, _Pr _Pred, const _Ty& _Val)
	{	
	struct _Unchecked_iterators { static void  _Deprecate(false_type) { } static void _Deprecate(true_type) { } }; (_Unchecked_iterators::_Deprecate(_Is_checked(_Dest)));
	return (_Replace_copy_if_no_deprecate(_First, _Last, _Dest, _Pred, _Val));
	}

 














		
template<class _FwdIt,
	class _Fn0> inline
	void _Generate_unchecked(_FwdIt _First, _FwdIt _Last, _Fn0& _Func)
	{	
	for (; _First != _Last; ++_First)
		*_First = _Func();
	}

template<class _FwdIt,
	class _Fn0> inline
	void generate(_FwdIt _First, _FwdIt _Last, _Fn0 _Func)
	{	
	;
	_Generate_unchecked(_Unchecked(_First), _Unchecked(_Last), _Func);
	}

		
template<class _OutIt,
	class _Diff,
	class _Fn0> inline
	_OutIt _Generate_n_unchecked(_OutIt _Dest, _Diff _Count, _Fn0& _Func)
	{	
	for (; 0 < _Count; --_Count, (void)++_Dest)
		*_Dest = _Func();
	return (_Dest);
	}

template<class _OutIt,
	class _Diff,
	class _Fn0> inline
	_OutIt generate_n(_OutIt _Dest, _Diff _Count, _Fn0 _Func)
	{	
	return (_Rechecked(_Dest,
		_Generate_n_unchecked(_Unchecked_n(_Dest, _Count), _Count, _Func)));
	}

 











		
template<class _InIt,
	class _OutIt,
	class _Ty> inline
	_OutIt _Remove_copy_unchecked(_InIt _First, _InIt _Last,
		_OutIt _Dest, const _Ty& _Val)
	{	
	for (; _First != _Last; ++_First)
		if (!(*_First == _Val))
			{	
			;
			*_Dest++ = *_First;
			}

	return (_Dest);
	}

template<class _InIt,
	class _OutIt,
	class _Ty> inline
	_OutIt remove_copy(_InIt _First, _InIt _Last,
		_OutIt _Dest, const _Ty& _Val)
	{	
	struct _Unchecked_iterators { static void  _Deprecate(false_type) { } static void _Deprecate(true_type) { } }; (_Unchecked_iterators::_Deprecate(_Is_checked(_Dest)));
	;
	return (_Rechecked(_Dest,
		_Remove_copy_unchecked(_Unchecked(_First), _Unchecked(_Last),
		_Unchecked_idl0(_Dest), _Val)));
	}

 













		
template<class _InIt,
	class _OutIt,
	class _Pr> inline
	_OutIt _Remove_copy_if_unchecked(_InIt _First, _InIt _Last,
		_OutIt _Dest, _Pr& _Pred)
	{	
	for (; _First != _Last; ++_First)
		if (!_Pred(*_First))
			{	
			;
			*_Dest++ = *_First;
			}

	return (_Dest);
	}

template<class _InIt,
	class _OutIt,
	class _Pr> inline
	_OutIt _Remove_copy_if_no_deprecate(_InIt _First, _InIt _Last,
		_OutIt _Dest, _Pr& _Pred)
	{	
	;
	return (_Rechecked(_Dest,
		_Remove_copy_if_unchecked(_Unchecked(_First), _Unchecked(_Last),
		_Unchecked_idl0(_Dest), _Pred)));
	}

template<class _InIt,
	class _OutIt,
	class _Pr> inline
	_OutIt remove_copy_if(_InIt _First, _InIt _Last,
		_OutIt _Dest, _Pr _Pred)
	{	
	struct _Unchecked_iterators { static void  _Deprecate(false_type) { } static void _Deprecate(true_type) { } }; (_Unchecked_iterators::_Deprecate(_Is_checked(_Dest)));
	return (_Remove_copy_if_no_deprecate(_First, _Last, _Dest, _Pred));
	}

 













		
template<class _FwdIt,
	class _Ty> inline
	_FwdIt _Remove_unchecked(_FwdIt _First, _FwdIt _Last, const _Ty& _Val)
	{	
	_First = _Find_unchecked(_First, _Last, _Val);
	_FwdIt _Next = _First;
	if (_First != _Last)
		{
		for (++_First; _First != _Last; ++_First)
			if (!(*_First == _Val))
				*_Next++ = ::std:: move(*_First);
		}

	return (_Next);
	}

template<class _FwdIt,
	class _Ty> inline
	_FwdIt remove(_FwdIt _First, _FwdIt _Last, const _Ty& _Val)
	{	
	;
	return (_Rechecked(_First,
		_Remove_unchecked(_Unchecked(_First), _Unchecked(_Last), _Val)));
	}

		
template<class _FwdIt,
	class _Pr> inline
	_FwdIt _Remove_if_unchecked(_FwdIt _First, _FwdIt _Last, _Pr& _Pred)
	{	
	_First = _Find_if_unchecked(_First, _Last, _Pred);
	_FwdIt _Next = _First;
	if (_First != _Last)
		{
		for (++_First; _First != _Last; ++_First)
		if (!_Pred(*_First))
			*_Next++ = ::std:: move(*_First);
		}

	return (_Next);
	}

template<class _FwdIt,
	class _Pr> inline
	_FwdIt remove_if(_FwdIt _First, _FwdIt _Last, _Pr _Pred)
	{	
	;
	return (_Rechecked(_First,
		_Remove_if_unchecked(_Unchecked(_First), _Unchecked(_Last), _Pred)));
	}

		
template<class _FwdIt,
	class _Pr> inline
	_FwdIt _Unique_unchecked(_FwdIt _First, _FwdIt _Last, _Pr& _Pred)
	{	
	if (_First != _Last)
		for (_FwdIt _Firstb; (void)(_Firstb = _First), ++_First != _Last; )
			if (_Pred(*_Firstb, *_First))
				{	
				for (; ++_First != _Last; )
					if (!_Pred(*_Firstb, *_First))
						*++_Firstb = ::std:: move(*_First);
				return (++_Firstb);
				}

	return (_Last);
	}

template<class _FwdIt,
	class _Pr> inline
	_FwdIt unique(_FwdIt _First, _FwdIt _Last, _Pr _Pred)
	{	
	;
	;
	return (_Rechecked(_First,
		_Unique_unchecked(_Unchecked(_First), _Unchecked(_Last), _Pred)));
	}

		
template<class _FwdIt> inline
	_FwdIt unique(_FwdIt _First, _FwdIt _Last)
	{	
	return (::std:: unique(_First, _Last, equal_to<>()));
	}

		
template<class _InIt,
	class _OutIt,
	class _Pr> inline
	_OutIt _Unique_copy_unchecked(_InIt _First, _InIt _Last,
		_OutIt _Dest, _Pr& _Pred, input_iterator_tag)
	{	
	if (_First != _Last)
		{
		_Iter_value_t<_InIt> _Val = *_First;

		for (*_Dest++ = _Val; ++_First != _Last; )
			if (!_Pred(_Val, *_First))
				{	
				_Val = *_First;
				*_Dest++ = _Val;
				}
		}

	return (_Dest);
	}

template<class _FwdIt,
	class _OutIt,
	class _Pr> inline
	_OutIt _Unique_copy_unchecked(_FwdIt _First, _FwdIt _Last,
		_OutIt _Dest, _Pr& _Pred, forward_iterator_tag)
	{	
	if (_First != _Last)
		{
		_FwdIt _Firstb = _First;

		for (*_Dest++ = *_Firstb; ++_First != _Last; )
			if (!_Pred(*_Firstb, *_First))
				{	
				_Firstb = _First;
				*_Dest++ = *_Firstb;
				}
		}

	return (_Dest);
	}

template<class _InIt,
	class _OutIt,
	class _Pr> inline
	_OutIt _Unique_copy_no_deprecate(_InIt _First, _InIt _Last,
		_OutIt _Dest, _Pr& _Pred)
	{	
	;
	;
	return (_Rechecked(_Dest,
		_Unique_copy_unchecked(_Unchecked(_First), _Unchecked(_Last),
			_Unchecked_idl0(_Dest), _Pred, _Iter_cat_t<_InIt>())));
	}

template<class _InIt,
	class _OutIt,
	class _Pr> inline
	_OutIt unique_copy(_InIt _First, _InIt _Last,
		_OutIt _Dest, _Pr _Pred)
	{	
	struct _Unchecked_iterators { static void  _Deprecate(false_type) { } static void _Deprecate(true_type) { } }; (_Unchecked_iterators::_Deprecate(_Is_checked(_Dest)));
	return (_Unique_copy_no_deprecate(_First, _Last, _Dest, _Pred));
	}

 













		
template<class _InIt,
	class _OutIt> inline
	_OutIt unique_copy(_InIt _First, _InIt _Last, _OutIt _Dest)
	{	
	return (::std:: unique_copy(_First, _Last, _Dest,
		equal_to<>()));
	}

 











		
template<class _BidIt,
	class _OutIt> inline
	_OutIt _Reverse_copy_unchecked(_BidIt _First, _BidIt _Last,
		_OutIt _Dest)
	{	
	for (; _First != _Last; ++_Dest)
		*_Dest = *--_Last;
	return (_Dest);
	}

template<class _BidIt,
	class _OutIt> inline
	_OutIt _Reverse_copy1(_BidIt _First, _BidIt _Last,
		_OutIt _Dest,
		bidirectional_iterator_tag, _Any_tag)
	{	
	return (_Rechecked(_Dest,
		_Reverse_copy_unchecked(_First, _Last, _Unchecked_idl0(_Dest))));
	}

template<class _BidIt,
	class _OutIt> inline
	_OutIt _Reverse_copy1(_BidIt _First, _BidIt _Last,
		_OutIt _Dest,
		random_access_iterator_tag, random_access_iterator_tag)
	{	
	;
	return (_Rechecked(_Dest,
		_Reverse_copy_unchecked(_First, _Last, _Unchecked(_Dest))));
	}

template<class _BidIt,
	class _OutIt> inline
	_OutIt reverse_copy(_BidIt _First, _BidIt _Last,
		_OutIt _Dest)
	{	
	struct _Unchecked_iterators { static void  _Deprecate(false_type) { } static void _Deprecate(true_type) { } }; (_Unchecked_iterators::_Deprecate(_Is_checked(_Dest)));
	;
	return (_Reverse_copy1(_Unchecked(_First), _Unchecked(_Last),
		_Dest, _Iter_cat_t<_BidIt>(), _Iter_cat_t<_OutIt>()));
	}

 












		
template<class _FwdIt,
	class _OutIt> inline
	_OutIt _Rotate_copy_unchecked(_FwdIt _First, _FwdIt _Mid, _FwdIt _Last,
		_OutIt _Dest)
	{	
	_Dest = _Copy_unchecked(_Mid, _Last, _Dest);
	return (_Copy_unchecked(_First, _Mid, _Dest));
	}

template<class _FwdIt,
	class _OutIt> inline
	_OutIt _Rotate_copy1(_FwdIt _First, _FwdIt _Mid, _FwdIt _Last,
		_OutIt _Dest, forward_iterator_tag, _Any_tag)
	{	
	return (_Rechecked(_Dest,
		_Rotate_copy_unchecked(_First, _Mid, _Last, _Unchecked_idl0(_Dest))));
	}

template<class _FwdIt,
	class _OutIt> inline
	_OutIt _Rotate_copy1(_FwdIt _First, _FwdIt _Mid, _FwdIt _Last,
		_OutIt _Dest, random_access_iterator_tag, random_access_iterator_tag)
	{	
	;
	return (_Rechecked(_Dest,
		_Rotate_copy_unchecked(_First, _Mid, _Last, _Unchecked(_Dest))));
	}

template<class _FwdIt,
	class _OutIt> inline
	_OutIt rotate_copy(_FwdIt _First, _FwdIt _Mid, _FwdIt _Last,
		_OutIt _Dest)
	{	
	;
	;
	return (_Rotate_copy1(_Unchecked(_First), _Unchecked(_Mid),
		_Unchecked(_Last), _Dest, _Iter_cat_t<_FwdIt>(), _Iter_cat_t<_OutIt>()));
	}

		
template<class _RanIt,
	class _Fn1> inline
	void _Random_shuffle_unchecked(_RanIt _First, _RanIt _Last, _Fn1& _Func)
	{	
	if (_Last - _First < 2)
		return;

	_RanIt _Next = _First;
	for (_Iter_diff_t<_RanIt> _Index = 2; ++_Next != _Last; ++_Index)
		{	
		_Iter_diff_t<_RanIt> _Off = _Func(_Index);

 










		::std:: iter_swap(_Next, _First + _Off);
		}
	}

template<class _RanIt,
	class _Fn1> inline
	void _Random_shuffle1(_RanIt _First, _RanIt _Last, _Fn1& _Func)
	{	
	;
	_Random_shuffle_unchecked(_Unchecked(_First), _Unchecked(_Last), _Func);
	}

template<class _RanIt,
	class _Urng> inline
	void shuffle(_RanIt _First, _RanIt _Last, _Urng&& _Func)
	{	
	typedef typename iterator_traits<_RanIt>::difference_type _Diff;
	typedef typename remove_reference<_Urng>::type _Urng0;
	_Rng_from_urng<_Diff, _Urng0> _Rng(_Func);
	_Random_shuffle1(_First, _Last, _Rng);
	}

 
		
template<class _RanIt,
	class _Fn1> inline
	void random_shuffle(_RanIt _First, _RanIt _Last, _Fn1&& _Func)
	{	
	_Random_shuffle1(_First, _Last, _Func);
	}

	
struct _Rand_urng_from_func
	{	
	typedef unsigned int result_type;

	static result_type (min)()
		{	
		return (0);
		}

	static result_type (max)()
		{	
		return (0x7fff);
		}

	result_type operator()()
		{	
		return (:: rand());
		}
	};

		
template<class _RanIt> inline
	void random_shuffle(_RanIt _First, _RanIt _Last)
	{	
	_Rand_urng_from_func _Func;
	::std:: shuffle(_First, _Last, _Func);
	}
 

		
template<class _FwdIt,
	class _Pr> inline
	_FwdIt _Partition_unchecked(_FwdIt _First, _FwdIt _Last, _Pr& _Pred,
		forward_iterator_tag)
	{	
	while (_First != _Last && _Pred(*_First))
		++_First;	

	if (_First == _Last)
		return (_First);	

	for (_FwdIt _Next = ::std:: next(_First); _Next != _Last; ++_Next)
		if (_Pred(*_Next))
			::std:: iter_swap(_First++, _Next);	

	return (_First);
	}

template<class _BidIt,
	class _Pr> inline
	_BidIt _Partition_unchecked(_BidIt _First, _BidIt _Last, _Pr& _Pred,
		bidirectional_iterator_tag)
	{	
	for (; ; ++_First)
		{	
		for (; _First != _Last && _Pred(*_First); ++_First)
			;	
		if (_First == _Last)
			break;	

		for (; _First != --_Last && !_Pred(*_Last); )
			;	
		if (_First == _Last)
			break;	

		::std:: iter_swap(_First, _Last);	
		}

	return (_First);
	}

template<class _FwdIt,
	class _Pr> inline
	_FwdIt partition(_FwdIt _First, _FwdIt _Last, _Pr _Pred)
	{	
	;
	return (_Rechecked(_First,
		_Partition_unchecked(_Unchecked(_First), _Unchecked(_Last), _Pred,
			_Iter_cat_t<_FwdIt>())));
	}

		
template<class _BidIt,
	class _Diff,
	class _Ty> inline
	_BidIt _Buffered_rotate_unchecked(_BidIt _First, _BidIt _Mid, _BidIt _Last,
		_Diff _Count1, _Diff _Count2, _Temp_iterator<_Ty>& _Tempbuf)
	{	
	if (_Count1 == 0 || _Count2 == 0)
		{	
		::std:: advance(_First, _Count2);
		return (_First);
		}
	else if (_Count1 <= _Count2 && _Count1 <= _Tempbuf._Maxlen())
		{	
		_Move_unchecked(_First, _Mid, _Tempbuf._Init());
		_Move_unchecked(_Mid, _Last, _First);
		return (_Move_backward_unchecked(_Tempbuf._First(), _Tempbuf._Last(),
			_Last));
		}
	else if (_Count2 <= _Tempbuf._Maxlen())
		{	
		_Move_unchecked(_Mid, _Last, _Tempbuf._Init());
		_Move_backward_unchecked(_First, _Mid, _Last);
		return (_Move_unchecked(_Tempbuf._First(), _Tempbuf._Last(), _First));
		}
	else
		{	
		return (_Rotate_unchecked(_First, _Mid, _Last));
		}
	}

template<class _BidIt,
	class _Pr,
	class _Diff,
	class _Ty> inline
	_BidIt _Stable_partition_unchecked1(_BidIt _First, _BidIt _Last, _Pr& _Pred,
		_Diff _Count, _Temp_iterator<_Ty>& _Tempbuf)
	{	
	if (_Count == 0)
		return (_First);
	else if (_Count == 1)
		return (_Pred(*_First) ? _Last : _First);
	else if (_Count <= _Tempbuf._Maxlen())
		{	
		_BidIt _Next = _First;
		for (_Tempbuf._Init(); _First != _Last; ++_First)
			if (_Pred(*_First))
				*_Next++ = ::std:: move(*_First);
			else
				*_Tempbuf++ = ::std:: move(*_First);

		_Move_unchecked(_Tempbuf._First(), _Tempbuf._Last(), _Next);	
		return (_Next);
		}
	else
		{	
		_BidIt _Mid = _First;
		::std:: advance(_Mid, _Count / 2);

		_BidIt _Left = _Stable_partition_unchecked1(_First, _Mid, _Pred,
			_Count / 2, _Tempbuf);	
		_BidIt _Right = _Stable_partition_unchecked1(_Mid, _Last, _Pred,
			_Count - _Count / 2, _Tempbuf);	

		_Diff _Count1 = ::std:: distance(_Left, _Mid);
		_Diff _Count2 = ::std:: distance(_Mid, _Right);

		return (_Buffered_rotate_unchecked(_Left, _Mid, _Right,
			_Count1, _Count2, _Tempbuf));	
		}
	}

template<class _BidIt,
	class _Pr> inline
	_BidIt _Stable_partition_unchecked(_BidIt _First, _BidIt _Last, _Pr& _Pred)
	{	
	if (_First == _Last)
		return (_First);
	_Iter_diff_t<_BidIt> _Count = ::std:: distance(_First, _Last);
	_Temp_iterator<_Iter_value_t<_BidIt>> _Tempbuf(_Count);
	return (_Stable_partition_unchecked1(_First, _Last, _Pred, _Count, _Tempbuf));
	}

template<class _BidIt,
	class _Pr> inline
	_BidIt stable_partition(_BidIt _First, _BidIt _Last, _Pr _Pred)
	{	
	;
	return (_Rechecked(_First,
		_Stable_partition_unchecked(_Unchecked(_First), _Unchecked(_Last), _Pred)));
	}

 





















  
 

		
template<class _RanIt,
	class _Diff,
	class _Ty,
	class _Pr> inline
	void _Push_heap_by_index(_RanIt _First, _Diff _Hole,
		_Diff _Top, _Ty&& _Val, _Pr& _Pred)
	{	
	for (_Diff _Idx = (_Hole - 1) / 2;
		_Top < _Hole && _Pred(*(_First + _Idx), _Val);
		_Idx = (_Hole - 1) / 2)
		{	
		*(_First + _Hole) = ::std:: move(*(_First + _Idx));
		_Hole = _Idx;
		}

	*(_First + _Hole) = ::std:: move(_Val);	
	}

template<class _RanIt,
	class _Pr> inline
	void _Push_heap_unchecked(_RanIt _First, _RanIt _Last, _Pr& _Pred)
	{	
	typedef _Iter_diff_t<_RanIt> _Diff;
	_Diff _Count = _Last - _First;
	if (2 <= _Count)
		{
		_Iter_value_t<_RanIt> _Val = ::std:: move(*--_Last);
		_Push_heap_by_index(_First, --_Count, _Diff(0), ::std:: move(_Val), _Pred);
		}
	}

template<class _RanIt,
	class _Pr> inline
	void push_heap(_RanIt _First, _RanIt _Last, _Pr _Pred)
	{	
	;
	_Push_heap_unchecked(_Unchecked(_First), _Unchecked(_Last), _Pred);
	}

		
template<class _RanIt> inline
	void push_heap(_RanIt _First, _RanIt _Last)
	{	
	::std:: push_heap(_First, _Last, less<>());
	}

		
template<class _RanIt,
	class _Diff,
	class _Ty,
	class _Pr> inline
	void _Pop_heap_hole_by_index(_RanIt _First, _Diff _Hole, _Diff _Bottom,
		_Ty&& _Val, _Pr& _Pred)
	{	
		
	const _Diff _Top = _Hole;
	_Diff _Idx = _Hole;

	
	
	const _Diff _Max_sequence_non_leaf = (_Bottom - 1) / 2;
	while (_Idx < _Max_sequence_non_leaf)
		{	
		_Idx = 2 * _Idx + 2;
		if (_Pred(*(_First + _Idx), *(_First + (_Idx - 1))))
			--_Idx;
		*(_First + _Hole) = ::std:: move(*(_First + _Idx));
		_Hole = _Idx;
		}

	if (_Idx == _Max_sequence_non_leaf && _Bottom % 2 == 0)
		{	
		*(_First + _Hole) = ::std:: move(*(_First + (_Bottom - 1)));
		_Hole = _Bottom - 1;
		}

	_Push_heap_by_index(_First, _Hole, _Top, ::std:: move(_Val), _Pred);
	}

template<class _RanIt,
	class _Ty,
	class _Pr> inline
	void _Pop_heap_hole_unchecked(_RanIt _First, _RanIt _Last, _RanIt _Dest,
		_Ty&& _Val, _Pr& _Pred)
	{	
		
		
	*_Dest = ::std:: move(*_First);
	_Pop_heap_hole_by_index(_First, _Iter_diff_t<_RanIt>(0), _Iter_diff_t<_RanIt>(_Last - _First),
		::std:: move(_Val), _Pred);
	}

template<class _RanIt,
	class _Pr> inline
	void _Pop_heap_unchecked(_RanIt _First, _RanIt _Last, _Pr& _Pred)
	{	
	if (2 <= _Last - _First)
		{
		--_Last;
		_Iter_value_t<_RanIt> _Val = ::std:: move(*_Last);
		_Pop_heap_hole_unchecked(_First, _Last, _Last,
			::std:: move(_Val), _Pred);
		}
	}

template<class _RanIt,
	class _Pr> inline
	void pop_heap(_RanIt _First, _RanIt _Last, _Pr _Pred)
	{	
	;
	;
	_Pop_heap_unchecked(_Unchecked(_First), _Unchecked(_Last), _Pred);
	}

		
template<class _RanIt> inline
	void pop_heap(_RanIt _First, _RanIt _Last)
	{	
	::std:: pop_heap(_First, _Last, less<>());
	}

		
template<class _RanIt,
	class _Pr> inline
	void _Make_heap_unchecked(_RanIt _First, _RanIt _Last, _Pr& _Pred)
	{	
	_Iter_diff_t<_RanIt> _Bottom = _Last - _First;
	if (2 <= _Bottom)
		{
		for (_Iter_diff_t<_RanIt> _Hole = _Bottom / 2; 0 < _Hole; )
			{	
			--_Hole;
			_Iter_value_t<_RanIt> _Val = ::std:: move(*(_First + _Hole));
			_Pop_heap_hole_by_index(_First, _Hole, _Bottom,
				::std:: move(_Val), _Pred);
			}
		}
	}

template<class _RanIt,
	class _Pr> inline
	void make_heap(_RanIt _First, _RanIt _Last, _Pr _Pred)
	{	
	;
	_Make_heap_unchecked(_Unchecked(_First), _Unchecked(_Last), _Pred);
	}

		
template<class _RanIt> inline
	void make_heap(_RanIt _First, _RanIt _Last)
	{	
	::std:: make_heap(_First, _Last, less<>());
	}

		
template<class _RanIt,
	class _Pr> inline
	void _Sort_heap_unchecked(_RanIt _First, _RanIt _Last, _Pr& _Pred)
	{	
	for (; 2 <= _Last - _First; --_Last)
		_Pop_heap_unchecked(_First, _Last, _Pred);
	}

template<class _RanIt,
	class _Pr> inline
	void sort_heap(_RanIt _First, _RanIt _Last, _Pr _Pred)
	{	
	;
	;
	;
	_Sort_heap_unchecked(_Unchecked(_First), _Unchecked(_Last), _Pred);
	}

		
template<class _RanIt> inline
	void sort_heap(_RanIt _First, _RanIt _Last)
	{	
	::std:: sort_heap(_First, _Last, less<>());
	}

		
template<class _FwdIt,
	class _Ty,
	class _Pr> inline
	_FwdIt _Lower_bound_unchecked(_FwdIt _First, _FwdIt _Last,
		const _Ty& _Val, _Pr& _Pred)
	{	
	_Iter_diff_t<_FwdIt> _Count = ::std:: distance(_First, _Last);

	while (0 < _Count)
		{	
		_Iter_diff_t<_FwdIt> _Count2 = _Count / 2;
		_FwdIt _Mid = _First;
		::std:: advance(_Mid, _Count2);

		if (_Pred(*_Mid, _Val))
			{	
			_First = ++_Mid;
			_Count -= _Count2 + 1;
			}
		else
			_Count = _Count2;
		}

	return (_First);
	}

template<class _FwdIt,
	class _Ty,
	class _Pr> inline
	_FwdIt lower_bound(_FwdIt _First, _FwdIt _Last,
		const _Ty& _Val, _Pr _Pred)
	{	
	;
	return (_Rechecked(_First,
		_Lower_bound_unchecked(_Unchecked(_First), _Unchecked(_Last), _Val, _Pred)));
	}

		
template<class _FwdIt,
	class _Ty> inline
	_FwdIt lower_bound(_FwdIt _First, _FwdIt _Last, const _Ty& _Val)
	{	
	return (::std:: lower_bound(_First, _Last, _Val, less<>()));
	}

		
template<class _FwdIt,
	class _Ty,
	class _Pr> inline
	_FwdIt _Upper_bound_unchecked(_FwdIt _First, _FwdIt _Last,
		const _Ty& _Val, _Pr& _Pred)
	{	
	_Iter_diff_t<_FwdIt> _Count = ::std:: distance(_First, _Last);

	while (0 < _Count)
		{	
		_Iter_diff_t<_FwdIt> _Count2 = _Count / 2;
		_FwdIt _Mid = _First;
		::std:: advance(_Mid, _Count2);

		if (!_Pred(_Val, *_Mid))
			{	
			_First = ++_Mid;
			_Count -= _Count2 + 1;
			}
		else
			_Count = _Count2;
		}

	return (_First);
	}

template<class _FwdIt,
	class _Ty,
	class _Pr> inline
	_FwdIt upper_bound(_FwdIt _First, _FwdIt _Last,
		const _Ty& _Val, _Pr _Pred)
	{	
	;
	return (_Rechecked(_First,
		_Upper_bound_unchecked(_Unchecked(_First), _Unchecked(_Last), _Val, _Pred)));
	}

		
template<class _FwdIt,
	class _Ty> inline
	_FwdIt upper_bound(_FwdIt _First, _FwdIt _Last, const _Ty& _Val)
	{	
	return (::std:: upper_bound(_First, _Last, _Val, less<>()));
	}

		
template<class _FwdIt,
	class _Ty,
	class _Pr> inline
	pair<_FwdIt, _FwdIt>
		_Equal_range_unchecked(_FwdIt _First, _FwdIt _Last,
			const _Ty& _Val, _Pr& _Pred)
	{	
	_Iter_diff_t<_FwdIt> _Count = ::std:: distance(_First, _Last);

	while (0 < _Count)
		{	
		_Iter_diff_t<_FwdIt> _Count2 = _Count / 2;
		_FwdIt _Mid = _First;
		::std:: advance(_Mid, _Count2);

		if (_Pred(*_Mid, _Val))
			{	
			_First = ++_Mid;
			_Count -= _Count2 + 1;
			}
		else if (_Pred(_Val, *_Mid))
			_Count = _Count2;	
		else
			{	
			_FwdIt _First2 = _Lower_bound_unchecked(_First, _Mid, _Val, _Pred);
			::std:: advance(_First, _Count);
			_FwdIt _Last2 = _Upper_bound_unchecked(++_Mid, _First, _Val, _Pred);
			return (pair<_FwdIt, _FwdIt>(_First2, _Last2));
			}
		}

	return (pair<_FwdIt, _FwdIt>(_First, _First));	
	}

template<class _FwdIt,
	class _Ty,
	class _Pr> inline
	pair<_FwdIt, _FwdIt>
		equal_range(_FwdIt _First, _FwdIt _Last,
			const _Ty& _Val, _Pr _Pred)
	{	
	;
	return (_Rechecked_both(_First, _Last,
		_Equal_range_unchecked(_Unchecked(_First), _Unchecked(_Last), _Val, _Pred)));
	}

		
template<class _FwdIt,
	class _Ty> inline
	pair<_FwdIt, _FwdIt>
		equal_range(_FwdIt _First, _FwdIt _Last,
			const _Ty& _Val)
	{	
	return (::std:: equal_range(_First, _Last, _Val, less<>()));
	}

		
template<class _FwdIt,
	class _Ty,
	class _Pr> inline
	bool _Binary_search_unchecked(_FwdIt _First, _FwdIt _Last,
		const _Ty& _Val, _Pr& _Pred)
	{	
	_First = _Lower_bound_unchecked(_First, _Last, _Val, _Pred);
	return (_First != _Last && !_Pred(_Val, *_First));
	}

template<class _FwdIt,
	class _Ty,
	class _Pr> inline
	bool binary_search(_FwdIt _First, _FwdIt _Last,
		const _Ty& _Val, _Pr _Pred)
	{	
	;
	return (_Binary_search_unchecked(_Unchecked(_First), _Unchecked(_Last),
		_Val, _Pred));
	}

		
template<class _FwdIt,
	class _Ty> inline
	bool binary_search(_FwdIt _First, _FwdIt _Last, const _Ty& _Val)
	{	
	return (::std:: binary_search(_First, _Last, _Val, less<>()));
	}

		
template<class _InIt1,
	class _InIt2,
	class _OutIt,
	class _Pr> inline
	_OutIt _Merge_unchecked(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _InIt2 _Last2,
		_OutIt _Dest, _Pr& _Pred)
	{	
	if (_First1 != _Last1 && _First2 != _Last2)
		for (; ; )
			{	
			if (_Pred(*_First2, *_First1))
				{
				*_Dest++ = *_First2++;
				if (_First2 == _Last2)
					break;
				}
			else
				{
				*_Dest++ = *_First1++;
				if (_First1 == _Last1)
					break;
				}
			}

	_Dest = _Copy_unchecked(_First1, _Last1, _Dest);	
	return (_Copy_unchecked(_First2, _Last2, _Dest));
	}

template<class _InIt1,
	class _InIt2,
	class _OutIt,
	class _Pr> inline
	_OutIt _Merge_no_deprecate1(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _InIt2 _Last2,
		_OutIt _Dest, _Pr& _Pred, input_iterator_tag,
		input_iterator_tag, _Any_tag)
	{	
	return (_Rechecked(_Dest,
		_Merge_unchecked(_First1, _Last1, _First2, _Last2, _Unchecked_idl0(_Dest), _Pred)));
	}

template<class _InIt1,
	class _InIt2,
	class _OutIt,
	class _Pr> inline
	_OutIt _Merge_no_deprecate1(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _InIt2 _Last2,
		_OutIt _Dest, _Pr& _Pred, random_access_iterator_tag,
		random_access_iterator_tag, random_access_iterator_tag)
	{	
	;
	return (_Rechecked(_Dest,
		_Merge_unchecked(_First1, _Last1, _First2, _Last2, _Unchecked(_Dest), _Pred)));
	}

template<class _InIt1,
	class _InIt2,
	class _OutIt,
	class _Pr> inline
	_OutIt _Merge_no_deprecate(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _InIt2 _Last2,
		_OutIt _Dest, _Pr& _Pred)
	{	
	;
	;
	;
	return (_Merge_no_deprecate1(_Unchecked(_First1), _Unchecked(_Last1),
		_Unchecked(_First2), _Unchecked(_Last2),
		_Dest, _Pred,
		_Iter_cat_t<_InIt1>(), _Iter_cat_t<_InIt2>(), _Iter_cat_t<_OutIt>()));
	}

template<class _InIt1,
	class _InIt2,
	class _OutIt,
	class _Pr> inline
	_OutIt merge(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _InIt2 _Last2,
		_OutIt _Dest, _Pr _Pred)
	{	
	struct _Unchecked_iterators { static void  _Deprecate(false_type) { } static void _Deprecate(true_type) { } }; (_Unchecked_iterators::_Deprecate(_Is_checked(_Dest)));
	return (_Merge_no_deprecate(_First1, _Last1, _First2, _Last2, _Dest, _Pred));
	}

 
















		
template<class _InIt1,
	class _InIt2,
	class _OutIt> inline
	_OutIt merge(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _InIt2 _Last2,
		_OutIt _Dest)
	{	
	return (::std:: merge(_First1, _Last1, _First2, _Last2, _Dest,
		less<>()));
	}

 













		
template<class _BidIt1,
	class _BidIt2,
	class _BidIt3,
	class _Pr> inline
	_BidIt3 _Buffered_merge_backward_unchecked(_BidIt1 _First1, _BidIt1 _Last1,
		_BidIt2 _First2, _BidIt2 _Last2, _BidIt3 _Dest, _Pr& _Pred,
		bool _In_place = false)
	{	
	if (_First1 != _Last1 && _First2 != _Last2)
		for (; ; )
			{	
			if (_Pred(*--_Last2, *--_Last1))
				{
				*--_Dest = ::std:: move(*_Last1);
				++_Last2;
				if (_First1 == _Last1)
					break;
				}
			else
				{
				*--_Dest = ::std:: move(*_Last2);
				++_Last1;
				if (_First2 == _Last2)
					break;
				}
			}

	_Dest = _Move_backward_unchecked(_First2, _Last2, _Dest);	
	if (!_In_place)
		_Dest = _Move_backward_unchecked(_First1, _Last1, _Dest);
	return (_Dest);
	}

template<class _InIt1,
	class _InIt2,
	class _OutIt,
	class _Pr> inline
	_OutIt _Buffered_merge_unchecked(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _InIt2 _Last2,
		_OutIt _Dest, _Pr& _Pred, bool _In_place = false)
	{	
	if (_First1 != _Last1 && _First2 != _Last2)
		for (; ; )
			{	
			if (_Pred(*_First2, *_First1))
				{	
				*_Dest++ = ::std:: move(*_First2++);
				if (_First2 == _Last2)
					break;
				}
			else
				{	
				*_Dest++ = ::std:: move(*_First1++);
				if (_First1 == _Last1)
					break;
				}
			}

	_Dest = _Move_unchecked(_First1, _Last1, _Dest);	
	if (!_In_place)
		_Dest = _Move_unchecked(_First2, _Last2, _Dest);
	return (_Dest);
	}

template<class _BidIt,
	class _Diff,
	class _Ty,
	class _Pr> inline
	void _Buffered_merge_unchecked(_BidIt _First, _BidIt _Mid, _BidIt _Last,
		_Diff _Count1, _Diff _Count2,
			_Temp_iterator<_Ty>& _Tempbuf, _Pr& _Pred)
	{	
	if (_Count1 == 0 || _Count2 == 0)
		;	
	else if (_Count1 + _Count2 == 2)
		{	
		if (_Pred(*_Mid, *_First))
			::std:: iter_swap(_First, _Mid);
		}
	else if (_Count1 <= _Count2 && _Count1 <= _Tempbuf._Maxlen())
		{	
		_Move_unchecked(_First, _Mid, _Tempbuf._Init());
		_Buffered_merge_unchecked(_Tempbuf._First(), _Tempbuf._Last(),
			_Mid, _Last, _First, _Pred, true);
		}
	else if (_Count2 <= _Tempbuf._Maxlen())
		{	
		_Move_unchecked(_Mid, _Last, _Tempbuf._Init());
		_Buffered_merge_backward_unchecked(_First, _Mid,
			_Tempbuf._First(), _Tempbuf._Last(), _Last, _Pred, true);
		}
	else
		{	
		_BidIt _Firstn, _Lastn;
		_Diff _Count1n, _Count2n;
		if (_Count2 < _Count1)
			{	
			_Count1n = _Count1 / 2;
			_Firstn = _First;
			::std:: advance(_Firstn, _Count1n);
			_Lastn = _Lower_bound_unchecked(_Mid, _Last, *_Firstn, _Pred);
			_Count2n = ::std:: distance(_Mid, _Lastn);
			}
		else
			{	
			_Count2n = _Count2 / 2;
			_Lastn = _Mid;
			::std:: advance(_Lastn, _Count2n);
			_Firstn = _Upper_bound_unchecked(_First, _Mid, *_Lastn, _Pred);
			_Count1n = ::std:: distance(_First, _Firstn);
			}
		_BidIt _Midn = _Buffered_rotate_unchecked(_Firstn, _Mid, _Lastn,
			_Count1 - _Count1n, _Count2n, _Tempbuf);	
		_Buffered_merge_unchecked(_First, _Firstn, _Midn,
			_Count1n, _Count2n, _Tempbuf, _Pred);	
		_Buffered_merge_unchecked(_Midn, _Lastn, _Last,
			_Count1 - _Count1n, _Count2 - _Count2n, _Tempbuf, _Pred);
		}
	}

template<class _BidIt,
	class _Pr> inline
	void _Inplace_merge_unchecked(_BidIt _First, _BidIt _Mid, _BidIt _Last, _Pr& _Pred)
	{	
	if (_First != _Mid && _Mid != _Last)
		{
		_Iter_diff_t<_BidIt> _Count1 = ::std:: distance(_First, _Mid);
		_Iter_diff_t<_BidIt> _Count2 = ::std:: distance(_Mid, _Last);
		_Temp_iterator<_Iter_value_t<_BidIt>> _Tempbuf(_Count1 < _Count2 ? _Count1 : _Count2);
		_Buffered_merge_unchecked(_First, _Mid, _Last,
			_Count1, _Count2, _Tempbuf, _Pred);
		}
	}

template<class _BidIt,
	class _Pr> inline
	void inplace_merge(_BidIt _First, _BidIt _Mid, _BidIt _Last, _Pr _Pred)
	{	
	;
	;
	_Inplace_merge_unchecked(
		_Unchecked(_First), _Unchecked(_Mid), _Unchecked(_Last), _Pred);
	}

		
template<class _BidIt> inline
	void inplace_merge(_BidIt _First, _BidIt _Mid, _BidIt _Last)
	{	
	::std:: inplace_merge(_First, _Mid, _Last, less<>());
	}

		
template<class _BidIt,
	class _Pr> inline
	void _Insertion_sort_unchecked(_BidIt _First, _BidIt _Last, _Pr& _Pred)
	{	
	if (_First != _Last)
		for (_BidIt _Next = _First; ++_Next != _Last; )
			{	
			_BidIt _Next1 = _Next;
			_Iter_value_t<_BidIt> _Val = ::std:: move(*_Next);

			if (_Pred(_Val, *_First))
				{	
				_Move_backward_unchecked(_First, _Next, ++_Next1);
				*_First = ::std:: move(_Val);
				}
			else
				{	
				for (_BidIt _First1 = _Next1;
					_Pred(_Val, *--_First1);
					_Next1 = _First1)
					*_Next1 = ::std:: move(*_First1);	
				*_Next1 = ::std:: move(_Val);	
				}
			}
	}

template<class _RanIt,
	class _Pr> inline
	void _Med3_unchecked(_RanIt _First, _RanIt _Mid, _RanIt _Last, _Pr& _Pred)
	{	
	if (_Pred(*_Mid, *_First))
		::std:: iter_swap(_Mid, _First);
	if (_Pred(*_Last, *_Mid))
		{	
		::std:: iter_swap(_Last, _Mid);
		if (_Pred(*_Mid, *_First))
			::std:: iter_swap(_Mid, _First);
		}
	}

template<class _RanIt,
	class _Pr> inline
	void _Guess_median_unchecked(_RanIt _First, _RanIt _Mid, _RanIt _Last, _Pr& _Pred)
	{	
	if (40 < _Last - _First)
		{	
		size_t _Step = (_Last - _First + 1) / 8;
		_Med3_unchecked(_First, _First + _Step, _First + 2 * _Step, _Pred);
		_Med3_unchecked(_Mid - _Step, _Mid, _Mid + _Step, _Pred);
		_Med3_unchecked(_Last - 2 * _Step, _Last - _Step, _Last, _Pred);
		_Med3_unchecked(_First + _Step, _Mid, _Last - _Step, _Pred);
		}
	else
		_Med3_unchecked(_First, _Mid, _Last, _Pred);
	}

template<class _RanIt,
	class _Pr> inline
	pair<_RanIt, _RanIt>
		_Partition_by_median_guess_unchecked(_RanIt _First, _RanIt _Last, _Pr& _Pred)
	{	
	_RanIt _Mid = _First + (_Last - _First) / 2;
	_Guess_median_unchecked(_First, _Mid, _Last - 1, _Pred);
	_RanIt _Pfirst = _Mid;
	_RanIt _Plast = _Pfirst + 1;

	while (_First < _Pfirst
		&& !_Pred(*(_Pfirst - 1), *_Pfirst)
		&& !_Pred(*_Pfirst, *(_Pfirst - 1)))
		--_Pfirst;
	while (_Plast < _Last
		&& !_Pred(*_Plast, *_Pfirst)
		&& !_Pred(*_Pfirst, *_Plast))
		++_Plast;

	_RanIt _Gfirst = _Plast;
	_RanIt _Glast = _Pfirst;

	for (; ; )
		{	
		for (; _Gfirst < _Last; ++_Gfirst)
			if (_Pred(*_Pfirst, *_Gfirst))
				;
			else if (_Pred(*_Gfirst, *_Pfirst))
				break;
			else if (_Plast++ != _Gfirst)
				::std:: iter_swap(_Plast - 1, _Gfirst);
		for (; _First < _Glast; --_Glast)
			if (_Pred(*(_Glast - 1), *_Pfirst))
				;
			else if (_Pred(*_Pfirst, *(_Glast - 1)))
				break;
			else if (--_Pfirst != _Glast - 1)
				::std:: iter_swap(_Pfirst, _Glast - 1);
		if (_Glast == _First && _Gfirst == _Last)
			return (pair<_RanIt, _RanIt>(_Pfirst, _Plast));

		if (_Glast == _First)
			{	
			if (_Plast != _Gfirst)
				::std:: iter_swap(_Pfirst, _Plast);
			++_Plast;
			::std:: iter_swap(_Pfirst++, _Gfirst++);
			}
		else if (_Gfirst == _Last)
			{	
			if (--_Glast != --_Pfirst)
				::std:: iter_swap(_Glast, _Pfirst);
			::std:: iter_swap(_Pfirst, --_Plast);
			}
		else
			::std:: iter_swap(_Gfirst++, --_Glast);
		}
	}

template<class _RanIt,
	class _Diff,
	class _Pr> inline
	void _Sort_unchecked1(_RanIt _First, _RanIt _Last, _Diff _Ideal, _Pr& _Pred)
	{	
	_Diff _Count;
	while (_ISORT_MAX < (_Count = _Last - _First) && 0 < _Ideal)
		{	
		pair<_RanIt, _RanIt> _Mid =
			_Partition_by_median_guess_unchecked(_First, _Last, _Pred);
		_Ideal /= 2, _Ideal += _Ideal / 2;	

		if (_Mid.first - _First < _Last - _Mid.second)
			{	
			_Sort_unchecked1(_First, _Mid.first, _Ideal, _Pred);
			_First = _Mid.second;
			}
		else
			{	
			_Sort_unchecked1(_Mid.second, _Last, _Ideal, _Pred);
			_Last = _Mid.first;
			}
		}

	if (_ISORT_MAX < _Count)
		{	
		_Make_heap_unchecked(_First, _Last, _Pred);
		_Sort_heap_unchecked(_First, _Last, _Pred);
		}
	else if (2 <= _Count)
		_Insertion_sort_unchecked(_First, _Last, _Pred);	
	}

template<class _RanIt,
	class _Pr> inline
	void _Sort_unchecked(_RanIt _First, _RanIt _Last, _Pr& _Pred)
	{	
	_Sort_unchecked1(_First, _Last, _Last - _First, _Pred);
	}

template<class _RanIt,
	class _Pr> inline
	void sort(_RanIt _First, _RanIt _Last, _Pr _Pred)
	{	
	;
	_Sort_unchecked(_Unchecked(_First), _Unchecked(_Last), _Pred);
	}

		
template<class _RanIt> inline
	void sort(_RanIt _First, _RanIt _Last)
	{	
	::std:: sort(_First, _Last, less<>());
	}

		
template<class _BidIt,
	class _OutIt,
	class _Diff,
	class _Pr> inline
	void _Chunked_merge_unchecked(_BidIt _First, _BidIt _Last, _OutIt _Dest,
		_Diff _Chunk, _Diff _Count, _Pr& _Pred)
	{	
	for (_Diff _Chunk2 = _Chunk * 2; _Chunk2 <= _Count; _Count -= _Chunk2)
		{	
		_BidIt _Mid1 = _First;
		::std:: advance(_Mid1, _Chunk);
		_BidIt _Mid2 = _Mid1;
		::std:: advance(_Mid2, _Chunk);

		_Dest = _Buffered_merge_unchecked(_First, _Mid1, _Mid1, _Mid2, _Dest, _Pred);
		_First = _Mid2;
		}

	if (_Count <= _Chunk)
		_Move_unchecked(_First, _Last, _Dest);	
	else
		{	
		_BidIt _Mid1 = _First;
		::std:: advance(_Mid1, _Chunk);

		_Buffered_merge_unchecked(_First, _Mid1, _Mid1, _Last, _Dest, _Pred);
		}
	}

template<class _BidIt,
	class _Diff,
	class _Ty,
	class _Pr> inline
	void _Buffered_merge_sort_unchecked(_BidIt _First, _BidIt _Last, _Diff _Count,
		_Temp_iterator<_Ty>& _Tempbuf, _Pr& _Pred)
	{	
	_BidIt _Mid = _First;
	for (_Diff _Nleft = _Count; _ISORT_MAX <= _Nleft; _Nleft -= _ISORT_MAX)
		{	
		_BidIt _Midn = _Mid;
		::std:: advance(_Midn, (int)_ISORT_MAX);

		_Insertion_sort_unchecked(_Mid, _Midn, _Pred);
		_Mid = _Midn;
		}
	_Insertion_sort_unchecked(_Mid, _Last, _Pred);	

	for (_Diff _Chunk = _ISORT_MAX; _Chunk < _Count; _Chunk *= 2)
		{	
		_Chunked_merge_unchecked(_First, _Last, _Tempbuf._Init(),
			_Chunk, _Count, _Pred);
		_Chunked_merge_unchecked(_Tempbuf._First(), _Tempbuf._Last(), _First,
			_Chunk *= 2, _Count, _Pred);
		}
	}

template<class _BidIt,
	class _Diff,
	class _Ty,
	class _Pr> inline
	void _Stable_sort_unchecked1(_BidIt _First, _BidIt _Last, _Diff _Count,
		_Temp_iterator<_Ty>& _Tempbuf, _Pr& _Pred)
	{	
	if (_Count <= _ISORT_MAX)
		_Insertion_sort_unchecked(_First, _Last, _Pred);	
	else
		{	
		_Diff _Count2 = (_Count + 1) / 2;
		_BidIt _Mid = _First;
		::std:: advance(_Mid, _Count2);

		if (_Count2 <= _Tempbuf._Maxlen())
			{	
			_Buffered_merge_sort_unchecked(_First, _Mid, _Count2, _Tempbuf, _Pred);
			_Buffered_merge_sort_unchecked(_Mid, _Last, _Count - _Count2,
				_Tempbuf, _Pred);
			}
		else
			{	
			_Stable_sort_unchecked1(_First, _Mid, _Count2, _Tempbuf, _Pred);
			_Stable_sort_unchecked1(_Mid, _Last, _Count - _Count2, _Tempbuf, _Pred);
			}

		_Buffered_merge_unchecked(_First, _Mid, _Last,
			_Count2, _Count - _Count2, _Tempbuf, _Pred);	
		}
	}

template<class _BidIt,
	class _Pr> inline
	void _Stable_sort_unchecked(_BidIt _First, _BidIt _Last, _Pr& _Pred)
	{	
	if (_First != _Last)
		{
		_Iter_diff_t<_BidIt> _Count = ::std:: distance(_First, _Last);
		_Temp_iterator<_Iter_value_t<_BidIt>> _Tempbuf((_Count + 1) / 2);
		_Stable_sort_unchecked1(_First, _Last, _Count, _Tempbuf, _Pred);
		}
	}

template<class _BidIt,
	class _Pr> inline
	void stable_sort(_BidIt _First, _BidIt _Last, _Pr _Pred)
	{	
	;
	_Stable_sort_unchecked(_Unchecked(_First), _Unchecked(_Last), _Pred);
	}

		
template<class _BidIt> inline
	void stable_sort(_BidIt _First, _BidIt _Last)
	{	
	::std:: stable_sort(_First, _Last, less<>());
	}

		
template<class _RanIt,
	class _Pr> inline
	void _Partial_sort_unchecked(_RanIt _First, _RanIt _Mid, _RanIt _Last,
		_Pr& _Pred)
	{	
	if (_First == _Mid)
		return;	
	_Make_heap_unchecked(_First, _Mid, _Pred);
	for (_RanIt _Next = _Mid; _Next < _Last; ++_Next)
		if (_Pred(*_Next, *_First))
			{	
			_Iter_value_t<_RanIt> _Val = ::std:: move(*_Next);
			_Pop_heap_hole_unchecked(_First, _Mid, _Next, ::std:: move(_Val), _Pred);
			}
	_Sort_heap_unchecked(_First, _Mid, _Pred);
	}

template<class _RanIt,
	class _Pr> inline
	void partial_sort(_RanIt _First, _RanIt _Mid, _RanIt _Last, _Pr _Pred)
	{	
	;
	;
	;
	_Partial_sort_unchecked(
		_Unchecked(_First), _Unchecked(_Mid), _Unchecked(_Last), _Pred);
	}

		
template<class _RanIt> inline
	void partial_sort(_RanIt _First, _RanIt _Mid, _RanIt _Last)
	{	
	::std:: partial_sort(_First, _Mid, _Last, less<>());
	}

		
template<class _InIt,
	class _RanIt,
	class _Pr> inline
	_RanIt _Partial_sort_copy_unchecked(_InIt _First1, _InIt _Last1,
		_RanIt _First2, _RanIt _Last2, _Pr& _Pred)
	{	
	_RanIt _Mid2 = _First2;
	if (_First1 != _Last1 && _First2 != _Last2)
		{
		for (; _First1 != _Last1 && _Mid2 != _Last2; ++_First1, (void)++_Mid2)
			*_Mid2 = *_First1;	
		_Make_heap_unchecked(_First2, _Mid2, _Pred);

		for (; _First1 != _Last1; ++_First1)
			if (_Pred(*_First1, *_First2))
				_Pop_heap_hole_by_index(_First2, _Iter_diff_t<_RanIt>(0), _Iter_diff_t<_RanIt>(_Mid2 - _First2),
					_Iter_value_t<_InIt>(*_First1), _Pred);	

		_Sort_heap_unchecked(_First2, _Mid2, _Pred);
		}

	return (_Mid2);
	}

template<class _InIt,
	class _RanIt,
	class _Pr> inline
	_RanIt partial_sort_copy(_InIt _First1, _InIt _Last1,
		_RanIt _First2, _RanIt _Last2, _Pr _Pred)
	{	
	;
	;
	return (_Rechecked(_First2,
		_Partial_sort_copy_unchecked(
			_Unchecked(_First1), _Unchecked(_Last1),
			_Unchecked(_First2), _Unchecked(_Last2), _Pred)));
	}

		
template<class _InIt,
	class _RanIt> inline
	_RanIt partial_sort_copy(_InIt _First1, _InIt _Last1,
		_RanIt _First2, _RanIt _Last2)
	{	
	return (::std:: partial_sort_copy(_First1, _Last1, _First2, _Last2,
		less<>()));
	}

		
template<class _RanIt,
	class _Pr> inline
	void _Nth_element_unchecked(_RanIt _First, _RanIt _Nth, _RanIt _Last, _Pr& _Pred)
	{	
	if (_Nth == _Last)
		return;	

	for (; _ISORT_MAX < _Last - _First; )
		{	
		pair<_RanIt, _RanIt> _Mid =
			_Partition_by_median_guess_unchecked(_First, _Last, _Pred);

		if (_Mid.second <= _Nth)
			_First = _Mid.second;
		else if (_Mid.first <= _Nth)
			return;	
		else
			_Last = _Mid.first;
		}

	_Insertion_sort_unchecked(_First, _Last, _Pred);	
	}

template<class _RanIt,
	class _Pr> inline
	void nth_element(_RanIt _First, _RanIt _Nth, _RanIt _Last, _Pr _Pred)
	{	
	;
	;
	;
	_Nth_element_unchecked(
		_Unchecked(_First), _Unchecked(_Nth), _Unchecked(_Last), _Pred);
	}

		
template<class _RanIt> inline
	void nth_element(_RanIt _First, _RanIt _Nth, _RanIt _Last)
	{	
	::std:: nth_element(_First, _Nth, _Last, less<>());
	}

		
template<class _InIt1,
	class _InIt2,
	class _Pr> inline
	bool _Includes_unchecked(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _InIt2 _Last2, _Pr& _Pred)
	{	
	for (; _First1 != _Last1 && _First2 != _Last2; )
		if (_Pred(*_First2, *_First1))
			return (false);
		else if (_Pred(*_First1, *_First2))
			++_First1;
		else
			{	
			++_First1;
			++_First2;
			}

	return (_First2 == _Last2);
	}

template<class _InIt1,
	class _InIt2,
	class _Pr> inline
	bool includes(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _InIt2 _Last2, _Pr _Pred)
	{	
	;
	;
	return (_Includes_unchecked(_Unchecked(_First1), _Unchecked(_Last1),
		_Unchecked(_First2), _Unchecked(_Last2), _Pred));
	}

		
template<class _InIt1,
	class _InIt2> inline
	bool includes(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _InIt2 _Last2)
	{	
	return (::std:: includes(_First1, _Last1, _First2, _Last2,
		less<>()));
	}

		
template<class _InIt1,
	class _InIt2,
	class _OutIt,
	class _Pr> inline
	_OutIt _Set_union_no_deprecate1(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _InIt2 _Last2, _OutIt _Dest, _Pr& _Pred)
	{	
	for (; _First1 != _Last1 && _First2 != _Last2; )
		if (_Pred(*_First1, *_First2))
			{	
			*_Dest++ = *_First1;
			++_First1;
			}
		else if (_Pred(*_First2, *_First1))
			{	
			*_Dest++ = *_First2;
			++_First2;
			}
		else
			{	
			*_Dest++ = *_First1;
			++_First1;
			++_First2;
			}
	_Dest = _Copy_no_deprecate(_First1, _Last1, _Dest);
	return (_Copy_no_deprecate(_First2, _Last2, _Dest));
	}

template<class _InIt1,
	class _InIt2,
	class _OutIt,
	class _Pr> inline
	_OutIt _Set_union_no_deprecate(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _InIt2 _Last2, _OutIt _Dest, _Pr& _Pred)
	{	
	;
	;
	;
	return (_Rechecked(_Dest,
		_Set_union_no_deprecate1(_Unchecked(_First1), _Unchecked(_Last1),
		_Unchecked(_First2), _Unchecked(_Last2),
		_Unchecked_idl0(_Dest), _Pred)));
	}

template<class _InIt1,
	class _InIt2,
	class _OutIt,
	class _Pr> inline
	_OutIt set_union(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _InIt2 _Last2, _OutIt _Dest, _Pr _Pred)
	{	
	struct _Unchecked_iterators { static void  _Deprecate(false_type) { } static void _Deprecate(true_type) { } }; (_Unchecked_iterators::_Deprecate(_Is_checked(_Dest)));
	return (_Set_union_no_deprecate(_First1, _Last1, _First2, _Last2, _Dest, _Pred));
	}

 















		
template<class _InIt1,
	class _InIt2,
	class _OutIt> inline
	_OutIt set_union(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _InIt2 _Last2, _OutIt _Dest)
	{	
	return (::std:: set_union(_First1, _Last1, _First2, _Last2, _Dest,
		less<>()));
	}

 












		
template<class _InIt1,
	class _InIt2,
	class _OutIt,
	class _Pr> inline
	_OutIt _Set_intersection_no_deprecate1(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _InIt2 _Last2, _OutIt _Dest, _Pr& _Pred)
	{	
	for (; _First1 != _Last1 && _First2 != _Last2; )
		if (_Pred(*_First1, *_First2))
			++_First1;
		else if (_Pred(*_First2, *_First1))
			++_First2;
		else
			{	
			;
			*_Dest++ = *_First1++;
			++_First2;
			}

	return (_Dest);
	}

template<class _InIt1,
	class _InIt2,
	class _OutIt,
	class _Pr> inline
	_OutIt _Set_intersection_no_deprecate(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _InIt2 _Last2, _OutIt _Dest, _Pr& _Pred)
	{	
	;
	;
	return (_Rechecked(_Dest,
		_Set_intersection_no_deprecate1(_Unchecked(_First1), _Unchecked(_Last1),
		_Unchecked(_First2), _Unchecked(_Last2),
		_Unchecked_idl0(_Dest), _Pred)));
	}

template<class _InIt1,
	class _InIt2,
	class _OutIt,
	class _Pr> inline
	_OutIt set_intersection(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _InIt2 _Last2, _OutIt _Dest, _Pr _Pred)
	{	
	struct _Unchecked_iterators { static void  _Deprecate(false_type) { } static void _Deprecate(true_type) { } }; (_Unchecked_iterators::_Deprecate(_Is_checked(_Dest)));
	return (_Set_intersection_no_deprecate(_First1, _Last1, _First2, _Last2, _Dest, _Pred));
	}

 















		
template<class _InIt1,
	class _InIt2,
	class _OutIt> inline
	_OutIt set_intersection(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _InIt2 _Last2, _OutIt _Dest)
	{	
	return (::std:: set_intersection(_First1, _Last1, _First2, _Last2, _Dest,
		less<>()));
	}

 












		
template<class _InIt1,
	class _InIt2,
	class _OutIt,
	class _Pr> inline
	_OutIt _Set_difference_no_deprecate1(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _InIt2 _Last2,
		_OutIt _Dest, _Pr& _Pred)
	{	
	for (; _First1 != _Last1 && _First2 != _Last2; )
		if (_Pred(*_First1, *_First2))
			{	
			;
			*_Dest++ = *_First1;
			++_First1;
			}
		else if (_Pred(*_First2, *_First1))
			++_First2;
		else
			{	
			++_First1;
			++_First2;
			}

	return (_Copy_no_deprecate(_First1, _Last1, _Dest));
	}

template<class _InIt1,
	class _InIt2,
	class _OutIt,
	class _Pr> inline
	_OutIt _Set_difference_no_deprecate(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _InIt2 _Last2,
		_OutIt _Dest, _Pr& _Pred)
	{	
	;
	;
	return (_Rechecked(_Dest,
		_Set_difference_no_deprecate1(_Unchecked(_First1), _Unchecked(_Last1),
		_Unchecked(_First2), _Unchecked(_Last2),
		_Unchecked_idl0(_Dest), _Pred)));
	}

template<class _InIt1,
	class _InIt2,
	class _OutIt,
	class _Pr> inline
	_OutIt set_difference(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _InIt2 _Last2,
		_OutIt _Dest, _Pr _Pred)
	{	
	struct _Unchecked_iterators { static void  _Deprecate(false_type) { } static void _Deprecate(true_type) { } }; (_Unchecked_iterators::_Deprecate(_Is_checked(_Dest)));
	return (_Set_difference_no_deprecate(_First1, _Last1, _First2, _Last2, _Dest, _Pred));
	}

 
















		
template<class _InIt1,
	class _InIt2,
	class _OutIt> inline
	_OutIt set_difference(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _InIt2 _Last2,
		_OutIt _Dest)
	{	
	return (::std:: set_difference(_First1, _Last1, _First2, _Last2, _Dest,
		less<>()));
	}

 













		
template<class _InIt1,
	class _InIt2,
	class _OutIt,
	class _Pr> inline
	_OutIt _Set_symmetric_difference_no_deprecate1(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _InIt2 _Last2,
		_OutIt _Dest, _Pr& _Pred)
	{	
	for (; _First1 != _Last1 && _First2 != _Last2; )
		if (_Pred(*_First1, *_First2))
			{	
			;
			*_Dest++ = *_First1;
			++_First1;
			}
		else if (_Pred(*_First2, *_First1))
			{	
			;
			*_Dest++ = *_First2;
			++_First2;
			}
		else
			{	
			++_First1;
			++_First2;
			}
	_Dest = _Copy_no_deprecate(_First1, _Last1, _Dest);
	return (_Copy_no_deprecate(_First2, _Last2, _Dest));
	}

template<class _InIt1,
	class _InIt2,
	class _OutIt,
	class _Pr> inline
	_OutIt _Set_symmetric_difference_no_deprecate(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _InIt2 _Last2,
		_OutIt _Dest, _Pr& _Pred)
	{	
	;
	;
	return (_Rechecked(_Dest,
		_Set_symmetric_difference_no_deprecate1(
		_Unchecked(_First1), _Unchecked(_Last1),
		_Unchecked(_First2), _Unchecked(_Last2),
		_Unchecked_idl0(_Dest), _Pred)));
	}

template<class _InIt1,
	class _InIt2,
	class _OutIt,
	class _Pr> inline
	_OutIt set_symmetric_difference(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _InIt2 _Last2,
		_OutIt _Dest, _Pr _Pred)
	{	
	struct _Unchecked_iterators { static void  _Deprecate(false_type) { } static void _Deprecate(true_type) { } }; (_Unchecked_iterators::_Deprecate(_Is_checked(_Dest)));
	return (_Set_symmetric_difference_no_deprecate(
		_First1, _Last1, _First2, _Last2, _Dest, _Pred));
	}

 
















		
template<class _InIt1,
	class _InIt2,
	class _OutIt> inline
	_OutIt set_symmetric_difference(_InIt1 _First1, _InIt1 _Last1,
		_InIt2 _First2, _InIt2 _Last2,
		_OutIt _Dest)
	{	
	return (::std:: set_symmetric_difference(_First1, _Last1, _First2, _Last2,
		_Dest, less<>()));
	}

 













		
template<class _FwdIt,
	class _Pr> inline
	_FwdIt _Max_element_unchecked(_FwdIt _First, _FwdIt _Last, _Pr& _Pred)
	{	
	_FwdIt _Found = _First;
	if (_First != _Last)
		for (; ++_First != _Last; )
			if (_Pred(*_Found, *_First))
				_Found = _First;
	return (_Found);
	}

template<class _FwdIt,
	class _Pr> inline
	_FwdIt max_element(_FwdIt _First, _FwdIt _Last, _Pr _Pred)
	{	
	;
	return (_Rechecked(_First,
		_Max_element_unchecked(_Unchecked(_First), _Unchecked(_Last), _Pred)));
	}

		
template<class _FwdIt> inline
	_FwdIt max_element(_FwdIt _First, _FwdIt _Last)
	{	
	return (::std:: max_element(_First, _Last, less<>()));
	}

		
template<class _FwdIt,
	class _Pr> inline
	_FwdIt _Min_element_unchecked(_FwdIt _First, _FwdIt _Last, _Pr& _Pred)
	{	
	_FwdIt _Found = _First;
	if (_First != _Last)
		for (; ++_First != _Last; )
			if (_Pred(*_First, *_Found))
				_Found = _First;
	return (_Found);
	}

template<class _FwdIt,
	class _Pr> inline
	_FwdIt min_element(_FwdIt _First, _FwdIt _Last, _Pr _Pred)
	{	
	;
	return (_Rechecked(_First,
		_Min_element_unchecked(_Unchecked(_First), _Unchecked(_Last), _Pred)));
	}

		
template<class _FwdIt> inline
	_FwdIt min_element(_FwdIt _First, _FwdIt _Last)
	{	
	return (::std:: min_element(_First, _Last, less<>()));
	}

		
template<class _FwdIt,
	class _Pr> inline
	pair<_FwdIt, _FwdIt>
		_Minmax_element_unchecked(_FwdIt _First, _FwdIt _Last, _Pr& _Pred)
	{	
	pair<_FwdIt, _FwdIt> _Found(_First, _First);

	if (_First != _Last)
		for (; ++_First != _Last; )
			{	
			_FwdIt _Next = _First;
			if (++_Next == _Last)
				{	
				if (_Pred(*_First, *_Found.first))
					_Found.first = _First;
				else if (!_Pred(*_First, *_Found.second))
					_Found.second = _First;
				}
			else
				{	
				if (_Pred(*_Next, *_First))
					{	
					if (_Pred(*_Next, *_Found.first))
						_Found.first = _Next;
					if (!_Pred(*_First, *_Found.second))
						_Found.second = _First;
					}
				else
					{	
					if (_Pred(*_First, *_Found.first))
						_Found.first = _First;
					if (!_Pred(*_Next, *_Found.second))
						_Found.second = _Next;
					}
				_First = _Next;
				}
			}

	return (_Found);
	}

template<class _FwdIt,
	class _Pr> inline
	pair<_FwdIt, _FwdIt>
		minmax_element(_FwdIt _First, _FwdIt _Last, _Pr _Pred)
	{	
	;
	return (_Rechecked_both(_First, _Last,
		_Minmax_element_unchecked(_Unchecked(_First), _Unchecked(_Last), _Pred)));
	}

		
template<class _FwdIt> inline
	pair<_FwdIt, _FwdIt>
		minmax_element(_FwdIt _First, _FwdIt _Last)
	{	
	return (::std:: minmax_element(_First, _Last, less<>()));
	}

		
template<class _Ty,
	class _Pr> inline
	constexpr const _Ty& (max)(const _Ty& _Left, const _Ty& _Right,
		_Pr _Pred)
		noexcept(noexcept(_Pred(_Left, _Right)))
	{	
	return (_Pred(_Left, _Right) ? _Right : _Left);
	}

template<class _Ty,
	class _Pr> inline
		
	_Ty (max)(::std:: initializer_list<_Ty> _Ilist, _Pr _Pred)
	{	
	const _Ty *_Res = _Max_element_unchecked(_Ilist.begin(), _Ilist.end(), _Pred);
	return (*_Res);
	}

		
template<class _Ty> inline

	 

	constexpr const _Ty& (max)(const _Ty& _Left, const _Ty& _Right)
		noexcept(noexcept(((_Left) < (_Right))))
	{	
	return (((_Left) < (_Right)) ? _Right : _Left);
	}

template<class _Ty> inline
		
	_Ty (max)(::std:: initializer_list<_Ty> _Ilist)
	{	
	return ((::std:: max)(_Ilist, less<>()));
	}

		
template<class _Ty,
	class _Pr> inline
	constexpr const _Ty& (min)(const _Ty& _Left, const _Ty& _Right,
		_Pr _Pred)
		noexcept(noexcept(_Pred(_Right, _Left)))
	{	
	return (_Pred(_Right, _Left) ? _Right : _Left);
	}

template<class _Ty,
	class _Pr> inline
		
	_Ty (min)(::std:: initializer_list<_Ty> _Ilist, _Pr _Pred)
	{	
	const _Ty *_Res = _Min_element_unchecked(_Ilist.begin(), _Ilist.end(), _Pred);
	return (*_Res);
	}

		
template<class _Ty> inline

	 

	constexpr const _Ty& (min)(const _Ty& _Left, const _Ty& _Right)
		noexcept(noexcept(((_Right) < (_Left))))
	{	
	return (((_Right) < (_Left)) ? _Right : _Left);
	}

template<class _Ty> inline
		
	_Ty (min)(::std:: initializer_list<_Ty> _Ilist)
	{	
	return ((::std:: min)(_Ilist, less<>()));
	}


		
template<class _Ty,
	class _Pr> inline
	constexpr pair<const _Ty&, const _Ty&>
		minmax(const _Ty& _Left, const _Ty& _Right, _Pr _Pred)
	{	
	return (_Pred(_Right, _Left)
		? pair<const _Ty&, const _Ty&>(_Right, _Left)
		: pair<const _Ty&, const _Ty&>(_Left, _Right));
	}

template<class _Ty,
	class _Pr> inline
		
	pair<_Ty, _Ty> minmax(::std:: initializer_list<_Ty> _Ilist,
		_Pr _Pred)
	{	
	pair<const _Ty *, const _Ty *> _Res = _Minmax_element_unchecked(
		_Ilist.begin(), _Ilist.end(), _Pred);
	return (pair<_Ty, _Ty>(*_Res.first, *_Res.second));
	}

		
template<class _Ty> inline
	constexpr pair<const _Ty&, const _Ty&>
		minmax(const _Ty& _Left, const _Ty& _Right)
	{	
	return (_Right < _Left
		? pair<const _Ty&, const _Ty&>(_Right, _Left)
		: pair<const _Ty&, const _Ty&>(_Left, _Right));
	}

template<class _Ty> inline
		
	pair<_Ty, _Ty> minmax(::std:: initializer_list<_Ty> _Ilist)
	{	
	return (::std:: minmax(_Ilist, less<>()));
	}

		
template<class _BidIt,
	class _Pr> inline
	bool _Next_permutation_unchecked(_BidIt _First, _BidIt _Last, _Pr& _Pred)
	{	
	_BidIt _Next = _Last;
	if (_First == _Last || _First == --_Next)
		return (false);

	for (; ; )
		{	
		_BidIt _Next1 = _Next;
		if (_Pred(*--_Next, *_Next1))
			{	
			_BidIt _Mid = _Last;
			for (; !_Pred(*_Next, *--_Mid); )
				;
			::std:: iter_swap(_Next, _Mid);
			_Reverse_unchecked(_Next1, _Last);
			return (true);
			}

		if (_Next == _First)
			{	
			_Reverse_unchecked(_First, _Last);
			return (false);
			}
		}
	}

template<class _BidIt,
	class _Pr> inline
	bool next_permutation(_BidIt _First, _BidIt _Last, _Pr _Pred)
	{	
	;
	return (_Next_permutation_unchecked(_Unchecked(_First), _Unchecked(_Last), _Pred));
	}

		
template<class _BidIt> inline
	bool next_permutation(_BidIt _First, _BidIt _Last)
	{	
	return (::std:: next_permutation(_First, _Last, less<>()));
	}

		
template<class _BidIt,
	class _Pr> inline
	bool _Prev_permutation_unchecked(_BidIt _First, _BidIt _Last, _Pr& _Pred)
	{	
	_BidIt _Next = _Last;
	if (_First == _Last || _First == --_Next)
		return (false);

	for (; ; )
		{	
		_BidIt _Next1 = _Next;
		if (_Pred(*_Next1, *--_Next))
			{	
			_BidIt _Mid = _Last;
			for (; !_Pred(*--_Mid, *_Next); )
				;
			::std:: iter_swap(_Next, _Mid);
			_Reverse_unchecked(_Next1, _Last);
			return (true);
			}

		if (_Next == _First)
			{	
			_Reverse_unchecked(_First, _Last);
			return (false);
			}
		}
	}

template<class _BidIt,
	class _Pr> inline
	bool prev_permutation(_BidIt _First, _BidIt _Last, _Pr _Pred)
	{	
	;
	return (_Prev_permutation_unchecked(_Unchecked(_First), _Unchecked(_Last), _Pred));
	}

		
template<class _BidIt> inline
	bool prev_permutation(_BidIt _First, _BidIt _Last)
	{	
	return (::std:: prev_permutation(_First, _Last, less<>()));
	}

		
template<class _RanIt,
	class _Pr> inline
	_RanIt _Is_heap_until_unchecked(_RanIt _First, _RanIt _Last, _Pr& _Pred)
	{	
	_Iter_diff_t<_RanIt> _Size = _Last - _First;

	if (2 <= _Size)
		for (_Iter_diff_t<_RanIt> _Off = 0; ++_Off < _Size; )
			if (_Pred(*(_First + (_Off - 1) / 2), *(_First + _Off)))
				return (_First + _Off);
	return (_Last);
	}

template<class _RanIt,
	class _Pr> inline
	_RanIt is_heap_until(_RanIt _First, _RanIt _Last, _Pr _Pred)
	{	
	;
	return (_Rechecked(_First,
		_Is_heap_until_unchecked(_Unchecked(_First), _Unchecked(_Last), _Pred)));
	}

template<class _RanIt,
	class _Pr> inline
	bool is_heap(_RanIt _First, _RanIt _Last, _Pr _Pred)
	{	
	;
	return (_Is_heap_until_unchecked(_Unchecked(_First), _Unchecked(_Last), _Pred) == _Unchecked(_Last));
	}

		
template<class _RanIt> inline
	_RanIt is_heap_until(_RanIt _First, _RanIt _Last)
	{	
	return (::std:: is_heap_until(_First, _Last, less<>()));
	}

template<class _RanIt> inline
	bool is_heap(_RanIt _First, _RanIt _Last)
	{	
	return (::std:: is_heap(_First, _Last, less<>()));
	}

		
template<class _FwdIt,
	class _Pr> inline
	_FwdIt _Is_sorted_until_unchecked(_FwdIt _First, _FwdIt _Last, _Pr& _Pred)
	{	
	if (_First != _Last)
		for (_FwdIt _Next = _First; ++_Next != _Last; ++_First)
			if (_Pred(*_Next, *_First))
				return (_Next);
	return (_Last);
	}

template<class _FwdIt,
	class _Pr> inline
	_FwdIt is_sorted_until(_FwdIt _First, _FwdIt _Last, _Pr _Pred)
	{	
	;
	return (_Rechecked(_First,
		_Is_sorted_until_unchecked(_Unchecked(_First), _Unchecked(_Last), _Pred)));
	}

template<class _FwdIt,
	class _Pr> inline
	bool is_sorted(_FwdIt _First, _FwdIt _Last, _Pr _Pred)
	{	
	;
	return (_Is_sorted_until_unchecked(_Unchecked(_First), _Unchecked(_Last), _Pred) == _Unchecked(_Last));
	}

		
template<class _FwdIt> inline
	_FwdIt is_sorted_until(_FwdIt _First, _FwdIt _Last)
	{	
	return (::std:: is_sorted_until(_First, _Last, less<>()));
	}

template<class _FwdIt> inline
	bool is_sorted(_FwdIt _First, _FwdIt _Last)
	{	
	return (::std:: is_sorted(_First, _Last, less<>()));
	}

































}
 
 #pragma warning(pop)
 #pragma pack(pop)











#pragma once





 #pragma pack(push,8)
 #pragma warning(push,3)
 
 
namespace std {
		
template<class _Kty,	
	class _Pr,	
	class _Alloc,	
	bool _Mfl>	
	class _Tset_traits
	{	
public:
	typedef _Kty key_type;
	typedef _Kty value_type;
	typedef _Pr key_compare;
	typedef _Alloc allocator_type;

	enum
		{	
		_Multi = _Mfl};

	typedef key_compare value_compare;

	static const _Kty& _Kfn(const value_type& _Val)
		{	
		return (_Val);
		}
	};

		
template<class _Kty,
	class _Pr = less<_Kty>,
	class _Alloc = allocator<_Kty> >
	class set
		: public _Tree<_Tset_traits<_Kty, _Pr, _Alloc, false> >
	{	
public:
	typedef set<_Kty, _Pr, _Alloc> _Myt;
	typedef _Tree<_Tset_traits<_Kty, _Pr, _Alloc, false> > _Mybase;
	typedef _Kty key_type;
	typedef _Pr key_compare;
	typedef typename _Mybase::value_compare value_compare;
	typedef typename _Mybase::allocator_type allocator_type;
	typedef typename _Mybase::size_type size_type;
	typedef typename _Mybase::difference_type difference_type;
	typedef typename _Mybase::pointer pointer;
	typedef typename _Mybase::const_pointer const_pointer;
	typedef typename _Mybase::reference reference;
	typedef typename _Mybase::const_reference const_reference;
	typedef typename _Mybase::iterator iterator;
	typedef typename _Mybase::const_iterator const_iterator;
	typedef typename _Mybase::reverse_iterator reverse_iterator;
	typedef typename _Mybase::const_reverse_iterator
		const_reverse_iterator;
	typedef typename _Mybase::value_type value_type;

	typedef typename _Mybase::_Alty _Alty;

	set()
		: _Mybase(key_compare())
		{	
		}

	explicit set(const allocator_type& _Al)
		: _Mybase(key_compare(), _Al)
		{	
		}

	set(const _Myt& _Right)
		: _Mybase(_Right,
			_Right._Getal().select_on_container_copy_construction())
		{	
		}

	set(const _Myt& _Right, const allocator_type& _Al)
		: _Mybase(_Right, _Al)
		{	
		}

	explicit set(const key_compare& _Pred)
		: _Mybase(_Pred)
		{	
		}

	set(const key_compare& _Pred, const allocator_type& _Al)
		: _Mybase(_Pred, _Al)
		{	
		}

	template<class _Iter>
		set(_Iter _First, _Iter _Last)
		: _Mybase(key_compare())
		{	
		this->insert(_First, _Last);
		}

	template<class _Iter>
		set(_Iter _First, _Iter _Last,
			const key_compare& _Pred)
		: _Mybase(_Pred)
		{	
		this->insert(_First, _Last);
		}

	template<class _Iter>
		set(_Iter _First, _Iter _Last,
			const key_compare& _Pred, const allocator_type& _Al)
		: _Mybase(_Pred, _Al)
		{	
		this->insert(_First, _Last);
		}

	_Myt& operator=(const _Myt& _Right)
		{	
		_Mybase::operator=(_Right);
		return (*this);
		}

	set(_Myt&& _Right)
		: _Mybase(::std:: move(_Right))
		{	
		}

	set(_Myt&& _Right, const allocator_type& _Al)
		: _Mybase(::std:: move(_Right), _Al)
		{	
		}

	_Myt& operator=(_Myt&& _Right)
		noexcept(_Alty::is_always_equal::value && is_nothrow_move_assignable<_Pr>::value)
		{	
		_Mybase::operator=(::std:: move(_Right));
		return (*this);
		}

	void swap(_Myt& _Right)
		noexcept(_Alty::is_always_equal::value && _Is_nothrow_swappable<_Pr>::value)
		{	
		_Mybase::swap(_Right);
		}

	set(::std:: initializer_list<value_type> _Ilist)
		: _Mybase(key_compare())
		{	
		this->insert(_Ilist);
		}

	set(::std:: initializer_list<value_type> _Ilist,
			const key_compare& _Pred)
		: _Mybase(_Pred)
		{	
		this->insert(_Ilist);
		}

	set(::std:: initializer_list<value_type> _Ilist,
			const key_compare& _Pred, const allocator_type& _Al)
		: _Mybase(_Pred, _Al)
		{	
		this->insert(_Ilist);
		}

	_Myt& operator=(::std:: initializer_list<value_type> _Ilist)
		{	
		this->clear();
		this->insert(_Ilist);
		return (*this);
		}
	};

template<class _Kty,
	class _Pr,
	class _Alloc> inline
	void swap(set<_Kty, _Pr, _Alloc>& _Left,
		set<_Kty, _Pr, _Alloc>& _Right)
		noexcept(noexcept(_Left.swap(_Right)))
	{	
	_Left.swap(_Right);
	}

		
template<class _Kty,
	class _Pr = less<_Kty>,
	class _Alloc = allocator<_Kty> >
	class multiset
		: public _Tree<_Tset_traits<_Kty, _Pr, _Alloc, true> >
	{	
public:
	typedef multiset<_Kty, _Pr, _Alloc> _Myt;
	typedef _Tree<_Tset_traits<_Kty, _Pr, _Alloc, true> > _Mybase;
	typedef _Kty key_type;
	typedef _Pr key_compare;
	typedef typename _Mybase::value_compare value_compare;
	typedef typename _Mybase::allocator_type allocator_type;
	typedef typename _Mybase::size_type size_type;
	typedef typename _Mybase::difference_type difference_type;
	typedef typename _Mybase::pointer pointer;
	typedef typename _Mybase::const_pointer const_pointer;
	typedef typename _Mybase::reference reference;
	typedef typename _Mybase::const_reference const_reference;
	typedef typename _Mybase::iterator iterator;
	typedef typename _Mybase::const_iterator const_iterator;
	typedef typename _Mybase::reverse_iterator reverse_iterator;
	typedef typename _Mybase::const_reverse_iterator
		const_reverse_iterator;
	typedef typename _Mybase::value_type value_type;

	typedef typename _Mybase::_Alty _Alty;

	multiset()
		: _Mybase(key_compare())
		{	
		}

	explicit multiset(const allocator_type& _Al)
		: _Mybase(key_compare(), _Al)
		{	
		}

	multiset(const _Myt& _Right)
		: _Mybase(_Right,
			_Right._Getal().select_on_container_copy_construction())
		{	
		}

	multiset(const _Myt& _Right, const allocator_type& _Al)
		: _Mybase(_Right, _Al)
		{	
		}

	explicit multiset(const key_compare& _Pred)
		: _Mybase(_Pred)
		{	
		}

	multiset(const key_compare& _Pred, const allocator_type& _Al)
		: _Mybase(_Pred, _Al)
		{	
		}

	template<class _Iter>
		multiset(_Iter _First, _Iter _Last)
		: _Mybase(key_compare())
		{	
		this->insert(_First, _Last);
		}

	template<class _Iter>
		multiset(_Iter _First, _Iter _Last,
			const key_compare& _Pred)
		: _Mybase(_Pred)
		{	
		this->insert(_First, _Last);
		}

	template<class _Iter>
		multiset(_Iter _First, _Iter _Last,
			const key_compare& _Pred, const allocator_type& _Al)
		: _Mybase(_Pred, _Al)
		{	
		this->insert(_First, _Last);
		}

	_Myt& operator=(const _Myt& _Right)
		{	
		_Mybase::operator=(_Right);
		return (*this);
		}

	multiset(_Myt&& _Right)
		: _Mybase(::std:: move(_Right))
		{	
		}

	multiset(_Myt&& _Right, const allocator_type& _Al)
		: _Mybase(::std:: move(_Right), _Al)
		{	
		}

	_Myt& operator=(_Myt&& _Right)
		noexcept(_Alty::is_always_equal::value && is_nothrow_move_assignable<_Pr>::value)
		{	
		_Mybase::operator=(::std:: move(_Right));
		return (*this);
		}

	template<class... _Valty>
		iterator emplace(_Valty&&... _Val)
		{	
		return (_Mybase::emplace(::std:: forward<_Valty>(_Val)...).first);
		}

	void swap(_Myt& _Right)
		noexcept(_Alty::is_always_equal::value && _Is_nothrow_swappable<_Pr>::value)
		{	
		_Mybase::swap(_Right);
		}

	multiset(::std:: initializer_list<value_type> _Ilist)
		: _Mybase(key_compare())
		{	
		this->insert(_Ilist);
		}

	multiset(::std:: initializer_list<value_type> _Ilist,
			const key_compare& _Pred)
		: _Mybase(_Pred)
		{	
		this->insert(_Ilist);
		}

	multiset(::std:: initializer_list<value_type> _Ilist,
			const key_compare& _Pred, const allocator_type& _Al)
		: _Mybase(_Pred, _Al)
		{	
		this->insert(_Ilist);
		}

	_Myt& operator=(::std:: initializer_list<value_type> _Ilist)
		{	
		this->clear();
		this->insert(_Ilist);
		return (*this);
		}
	};

template<class _Kty,
	class _Pr,
	class _Alloc> inline
	void swap(multiset<_Kty, _Pr, _Alloc>& _Left,
		multiset<_Kty, _Pr, _Alloc>& _Right)
		noexcept(noexcept(_Left.swap(_Right)))
	{	
	_Left.swap(_Right);
	}
}
 
 #pragma warning(pop)
 #pragma pack(pop)





































































































































































































































#pragma once












































































































































 





























































#pragma warning(push)

#pragma warning(disable:4001) 


#pragma once





























































































































































































































































#pragma warning(pop)













#pragma once


#pragma region Application Family





















































































































#pragma warning(disable:4514)

#pragma warning(disable:4103)


#pragma warning(push)

#pragma warning(disable:4001)
#pragma warning(disable:4201)
#pragma warning(disable:4214)










#pragma once




































































































































































































































































































































__pragma(pack(push, 8)) extern "C" {




typedef enum _EXCEPTION_DISPOSITION
{
    ExceptionContinueExecution,
    ExceptionContinueSearch,
    ExceptionNestedException,
    ExceptionCollidedUnwind
} EXCEPTION_DISPOSITION;

















    

        struct _EXCEPTION_RECORD;
        struct _CONTEXT;
        struct _DISPATCHER_CONTEXT;

        __declspec(dllimport) EXCEPTION_DISPOSITION __C_specific_handler(
                 struct _EXCEPTION_RECORD*   ExceptionRecord,
                 void*                       EstablisherFrame,
              struct _CONTEXT*            ContextRecord,
              struct _DISPATCHER_CONTEXT* DispatcherContext
            );

    












unsigned long __cdecl _exception_code(void);
void *        __cdecl _exception_info(void);
int           __cdecl _abnormal_termination(void);










} __pragma(pack(pop))









#pragma once




































































































































































































































































































































__pragma(pack(push, 8)) extern "C" {










} __pragma(pack(pop))

















#pragma once


extern "C" {


















#pragma once


















#pragma once




 


  
 

 



  
 

 
  
  
 
























extern "C" {


































































































































































































































































































































































































































                                

}























































































































































































































#pragma once






























































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































#pragma once



extern "C" {







































































































    
    
    
    
    
    
    
    
    
    
    




    
    
    
    
    

    
    
    
    
    
    
    

    
    
    
    



    
    


    
    
    
    
    
    
    
    
    
    
    
    


    
    


    
    


    
    



    
    









    
    
    
    
    
    
    
    
    
    

    
    
    
    
    
    
    
    
    
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    

    
    


    


    
    
    
    
    

    


    
    
    
    
    

    


    
    
    
    
    

    


    
    
    
    
    


    




    
    
    
    
    

    


    
    
    
    
    


    


    
    
    
    
    
    

    


    
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    
    
    
    

    
    
    
    
    

    
    
    
    
    

    
    
    
    
    

    
    
    
    
    

    
    
    
    
    

    
    
    
    
    
    

    
    
    
    
    

    
    
    
    
    
    

    
    
    
    
    
    

    
    
    
    
    

    
    
    
    
    
    

    
    
    
    
    
    

    
    

    
    
    
    

    

    
    
    

    

    
    
    
    
    
    
    
    
    
    
    
    

    
    


    
    
    
    
    
    
    
    
    
    
    

    
    
    
    
    

    
    
    
    
    
    


    
    
    
    
     
    

    
    
        
        
        
        
    
    
    
    
    
    

    
    
    

    
    
    
    
    
    
    
    
    
    
    
    

    
    
    
    
    
    

    

    
    
    
    
    
    
    
    
    
    
    
    
    


    
    
    
    
    
    


    
    
    
    
    
    
    
    
    
    
    
    
    
    

    
    
    
    
    

    
    
    
    
    
    

    
    
    
    
    
    
    
    
    

    
    

    
    

    
    
    
    
    
    


    
    
    

    
	

    
    


    
    


    
    
    


    
    
    


    
    











































































    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    



    
    
    

    
    
    


}




























#pragma region Application Family























extern "C" {








typedef unsigned long ULONG;
typedef ULONG *PULONG;
typedef unsigned short USHORT;
typedef USHORT *PUSHORT;
typedef unsigned char UCHAR;
typedef UCHAR *PUCHAR;
typedef   char *PSZ;


































































































typedef unsigned long       DWORD;
typedef int                 BOOL;
typedef unsigned char       BYTE;
typedef unsigned short      WORD;
typedef float               FLOAT;
typedef FLOAT               *PFLOAT;
typedef BOOL            *PBOOL;
typedef BOOL             *LPBOOL;
typedef BYTE            *PBYTE;
typedef BYTE             *LPBYTE;
typedef int             *PINT;
typedef int              *LPINT;
typedef WORD            *PWORD;
typedef WORD             *LPWORD;
typedef long             *LPLONG;
typedef DWORD           *PDWORD;
typedef DWORD            *LPDWORD;
typedef void             *LPVOID;
typedef const void       *LPCVOID;

typedef int                 INT;
typedef unsigned int        UINT;
typedef unsigned int        *PUINT;

























#pragma warning(push)

#pragma warning(disable:4201) 
#pragma warning(disable:4214) 


extern "C" {










#pragma once





__pragma(pack(push, 8)) extern "C" {









  __declspec(dllimport) int __cdecl _isctype(  int _C,   int _Type);
  __declspec(dllimport) int __cdecl _isctype_l(  int _C,   int _Type,   _locale_t _Locale);
   __declspec(dllimport) int __cdecl isalpha(  int _C);
  __declspec(dllimport) int __cdecl _isalpha_l(  int _C,   _locale_t _Locale);
   __declspec(dllimport) int __cdecl isupper(  int _C);
  __declspec(dllimport) int __cdecl _isupper_l(  int _C,   _locale_t _Locale);
   __declspec(dllimport) int __cdecl islower(  int _C);
  __declspec(dllimport) int __cdecl _islower_l(  int _C,   _locale_t _Locale);


   __declspec(dllimport) int __cdecl isdigit(  int _C);

  __declspec(dllimport) int __cdecl _isdigit_l(  int _C,   _locale_t _Locale);
  __declspec(dllimport) int __cdecl isxdigit(  int _C);
  __declspec(dllimport) int __cdecl _isxdigit_l(  int _C,   _locale_t _Locale);


   __declspec(dllimport) int __cdecl isspace(  int _C);

  __declspec(dllimport) int __cdecl _isspace_l(  int _C,   _locale_t _Locale);
  __declspec(dllimport) int __cdecl ispunct(  int _C);
  __declspec(dllimport) int __cdecl _ispunct_l(  int _C,   _locale_t _Locale);
  __declspec(dllimport) int __cdecl isblank(  int _C);
  __declspec(dllimport) int __cdecl _isblank_l(  int _C,   _locale_t _Locale);
   __declspec(dllimport) int __cdecl isalnum(  int _C);
  __declspec(dllimport) int __cdecl _isalnum_l(  int _C,   _locale_t _Locale);
  __declspec(dllimport) int __cdecl isprint(  int _C);
  __declspec(dllimport) int __cdecl _isprint_l(  int _C,   _locale_t _Locale);
  __declspec(dllimport) int __cdecl isgraph(  int _C);
  __declspec(dllimport) int __cdecl _isgraph_l(  int _C,   _locale_t _Locale);
  __declspec(dllimport) int __cdecl iscntrl(  int _C);
  __declspec(dllimport) int __cdecl _iscntrl_l(  int _C,   _locale_t _Locale);


   __declspec(dllimport) int __cdecl toupper(  int _C);


   __declspec(dllimport) int __cdecl tolower(  int _C);

   __declspec(dllimport) int __cdecl _tolower(  int _C);
  __declspec(dllimport) int __cdecl _tolower_l(  int _C,   _locale_t _Locale);
   __declspec(dllimport) int __cdecl _toupper(  int _C);
  __declspec(dllimport) int __cdecl _toupper_l(  int _C,   _locale_t _Locale);

  __declspec(dllimport) int __cdecl __isascii(  int _C);
  __declspec(dllimport) int __cdecl __toascii(  int _C);
  __declspec(dllimport) int __cdecl __iscsymf(  int _C);
  __declspec(dllimport) int __cdecl __iscsym(  int _C);









 
    
    
    






















    
    
    



        
    



    
    
    
    
    
    
    
    



    














    __inline __crt_locale_data_public* __cdecl __acrt_get_locale_data_prefix(void const volatile* const _LocalePointers)
    {
        _locale_t const _TypedLocalePointers = (_locale_t)_LocalePointers;
        return (__crt_locale_data_public*)_TypedLocalePointers->locinfo;
    }

    



    __inline int __cdecl _chvalidchk_l(
              int       const _C,
              int       const _Mask,
          _locale_t const _Locale
        )
    {
        


        if (_Locale)
        {
            return __acrt_get_locale_data_prefix(_Locale)->_locale_pctype[_C] & _Mask;
        }
            
        return (__pctype_func()[(_C)] & (_Mask));
        
    }

    
    

    __inline int __cdecl _ischartype_l(
              int       const _C,
              int       const _Mask,
          _locale_t const _Locale
        )
    {
        if (_Locale && __acrt_get_locale_data_prefix(_Locale)->_locale_mb_cur_max > 1)
        {
            return _isctype_l(_C, _Mask, _Locale);
        }

        return _chvalidchk_l(_C, _Mask, _Locale);
    }

    
    
    
    
    
    
    
    
    
    
    
    

    
    

    
    


    
    
    
    
    

    
    
    
    





    
    
    
    
    





} __pragma(pack(pop))



























































































#pragma once



extern "C" {


    
    
    

    
    
    

    


    
    
    
    


    



    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    
    
    
    
    
    
    


    
    
    
    
    
    

    
    

    
    
    
    
    
    

    
    

    
    
    
    
    
    

    
    

    
    
    
    
    
    

    
    

    
    
    
    
    
    

    
    

    
    
    
    
    
    
    

    
    

    
    
    
    
    
    

    
    

    
    
    
    
    
    
    

    
    

    
    
    
    
    
    
    

    
    

    
    
    
    
    
    

    
    

    
    
    
    
    
    
    

    
    

    
    
    
    
    
    
    

    
    















}


















































































































 
 typedef unsigned __int64 POINTER_64_INT;
 
  
 



































#pragma once



extern "C" {


typedef signed char         INT8, *PINT8;
typedef signed short        INT16, *PINT16;
typedef signed int          INT32, *PINT32;
typedef signed __int64      INT64, *PINT64;
typedef unsigned char       UINT8, *PUINT8;
typedef unsigned short      UINT16, *PUINT16;
typedef unsigned int        UINT32, *PUINT32;
typedef unsigned __int64    UINT64, *PUINT64;





typedef signed int LONG32, *PLONG32;





typedef unsigned int ULONG32, *PULONG32;
typedef unsigned int DWORD32, *PDWORD32;





























    typedef __int64 INT_PTR, *PINT_PTR;
    typedef unsigned __int64 UINT_PTR, *PUINT_PTR;

    typedef __int64 LONG_PTR, *PLONG_PTR;
    typedef unsigned __int64 ULONG_PTR, *PULONG_PTR;

    























typedef __int64 SHANDLE_PTR;
typedef unsigned __int64 HANDLE_PTR;
typedef unsigned int UHALF_PTR, *PUHALF_PTR;
typedef int HALF_PTR, *PHALF_PTR;


__inline
unsigned long
HandleToULong(
    const void *h
    )
{
    return((unsigned long) (ULONG_PTR) h );
}

__inline
long
HandleToLong(
    const void *h
    )
{
    return((long) (LONG_PTR) h );
}

__inline
void *
ULongToHandle(
    const unsigned long h
    )
{
    return((void *) (UINT_PTR) h );
}


__inline
void *
LongToHandle(
    const long h
    )
{
    return((void *) (INT_PTR) h );
}


__inline
unsigned long
PtrToUlong(
    const void  *p
    )
{
    return((unsigned long) (ULONG_PTR) p );
}

__inline
unsigned int
PtrToUint(
    const void  *p
    )
{
    return((unsigned int) (UINT_PTR) p );
}

__inline
unsigned short
PtrToUshort(
    const void  *p
    )
{
    return((unsigned short) (unsigned long) (ULONG_PTR) p );
}

__inline
long
PtrToLong(
    const void  *p
    )
{
    return((long) (LONG_PTR) p );
}

__inline
int
PtrToInt(
    const void  *p
    )
{
    return((int) (INT_PTR) p );
}

__inline
short
PtrToShort(
    const void  *p
    )
{
    return((short) (long) (LONG_PTR) p );
}

__inline
void *
IntToPtr(
    const int i
    )

{
    return( (void *)(INT_PTR)i );
}

__inline
void *
UIntToPtr(
    const unsigned int ui
    )

{
    return( (void *)(UINT_PTR)ui );
}

__inline
void *
LongToPtr(
    const long l
    )

{
    return( (void *)(LONG_PTR)l );
}

__inline
void *
ULongToPtr(
    const unsigned long ul
    )

{
    return( (void *)(ULONG_PTR)ul );
}






__inline
void *
Ptr32ToPtr(
    const void * __ptr32 p
    )
{
    return((void *) (ULONG_PTR) (unsigned long) p);
}

__inline
void *
Handle32ToHandle(
    const void * __ptr32 h
    )
{
    return((void *) (LONG_PTR) (long) h);
}

__inline
void * __ptr32
PtrToPtr32(
    const void *p
    )
{
    return((void * __ptr32) (unsigned long) (ULONG_PTR) p);
}
































































































typedef ULONG_PTR SIZE_T, *PSIZE_T;
typedef LONG_PTR SSIZE_T, *PSSIZE_T;















































typedef ULONG_PTR DWORD_PTR, *PDWORD_PTR;





typedef __int64 LONG64, *PLONG64;






typedef unsigned __int64 ULONG64, *PULONG64;
typedef unsigned __int64 DWORD64, *PDWORD64;







typedef ULONG_PTR KAFFINITY;
typedef KAFFINITY *PKAFFINITY;




}













































































































































































typedef void *PVOID;
typedef void * __ptr64 PVOID64;








































typedef char CHAR;
typedef short SHORT;
typedef long LONG;

typedef int INT;








typedef wchar_t WCHAR;    





typedef WCHAR *PWCHAR, *LPWCH, *PWCH;
typedef const WCHAR *LPCWCH, *PCWCH;

typedef   WCHAR *NWPSTR, *LPWSTR, *PWSTR;
typedef   PWSTR *PZPWSTR;
typedef   const PWSTR *PCZPWSTR;
typedef   WCHAR __unaligned *LPUWSTR, *PUWSTR;
typedef   const WCHAR *LPCWSTR, *PCWSTR;
typedef   PCWSTR *PZPCWSTR;
typedef   const PCWSTR *PCZPCWSTR;
typedef   const WCHAR __unaligned *LPCUWSTR, *PCUWSTR;

typedef   WCHAR *PZZWSTR;
typedef   const WCHAR *PCZZWSTR;
typedef   WCHAR __unaligned *PUZZWSTR;
typedef   const WCHAR __unaligned *PCUZZWSTR;

typedef  WCHAR *PNZWCH;
typedef  const WCHAR *PCNZWCH;
typedef  WCHAR __unaligned *PUNZWCH;
typedef  const WCHAR __unaligned *PCUNZWCH;



typedef const WCHAR *LPCWCHAR, *PCWCHAR;
typedef const WCHAR __unaligned *LPCUWCHAR, *PCUWCHAR;





typedef unsigned long UCSCHAR;



















typedef UCSCHAR *PUCSCHAR;
typedef const UCSCHAR *PCUCSCHAR;

typedef UCSCHAR *PUCSSTR;
typedef UCSCHAR __unaligned *PUUCSSTR;

typedef const UCSCHAR *PCUCSSTR;
typedef const UCSCHAR __unaligned *PCUUCSSTR;

typedef UCSCHAR __unaligned *PUUCSCHAR;
typedef const UCSCHAR __unaligned *PCUUCSCHAR;







typedef CHAR *PCHAR, *LPCH, *PCH;
typedef const CHAR *LPCCH, *PCCH;

typedef   CHAR *NPSTR, *LPSTR, *PSTR;
typedef   PSTR *PZPSTR;
typedef   const PSTR *PCZPSTR;
typedef   const CHAR *LPCSTR, *PCSTR;
typedef   PCSTR *PZPCSTR;
typedef   const PCSTR *PCZPCSTR;

typedef   CHAR *PZZSTR;
typedef   const CHAR *PCZZSTR;

typedef  CHAR *PNZCH;
typedef  const CHAR *PCNZCH;

































typedef char TCHAR, *PTCHAR;
typedef unsigned char TBYTE , *PTBYTE ;



typedef LPCH LPTCH, PTCH;
typedef LPCCH LPCTCH, PCTCH;
typedef LPSTR PTSTR, LPTSTR, PUTSTR, LPUTSTR;
typedef LPCSTR PCTSTR, LPCTSTR, PCUTSTR, LPCUTSTR;
typedef PZZSTR PZZTSTR, PUZZTSTR;
typedef PCZZSTR PCZZTSTR, PCUZZTSTR;
typedef PZPSTR PZPTSTR;
typedef PNZCH PNZTCH, PUNZTCH;
typedef PCNZCH PCNZTCH, PCUNZTCH;






typedef SHORT *PSHORT;  
typedef LONG *PLONG;    








typedef struct _PROCESSOR_NUMBER {
    WORD   Group;
    BYTE  Number;
    BYTE  Reserved;
} PROCESSOR_NUMBER, *PPROCESSOR_NUMBER;






typedef struct _GROUP_AFFINITY {
    KAFFINITY Mask;
    WORD   Group;
    WORD   Reserved[3];
} GROUP_AFFINITY, *PGROUP_AFFINITY;








typedef void *HANDLE;









typedef HANDLE *PHANDLE;







typedef BYTE   FCHAR;
typedef WORD   FSHORT;
typedef DWORD  FLONG;










typedef   long HRESULT;






    























































typedef char CCHAR;          
typedef DWORD LCID;         
typedef PDWORD PLCID;       
typedef WORD   LANGID;      








typedef enum {
    UNSPECIFIED_COMPARTMENT_ID = 0,
    DEFAULT_COMPARTMENT_ID
} COMPARTMENT_ID, *PCOMPARTMENT_ID;



























typedef struct _FLOAT128 {
    __int64 LowPart;
    __int64 HighPart;
} FLOAT128;

typedef FLOAT128 *PFLOAT128;









typedef __int64 LONGLONG;
typedef unsigned __int64 ULONGLONG;




















typedef LONGLONG *PLONGLONG;
typedef ULONGLONG *PULONGLONG;



typedef LONGLONG USN;




typedef union _LARGE_INTEGER {
    struct {
        DWORD LowPart;
        LONG HighPart;
    } ;
    struct {
        DWORD LowPart;
        LONG HighPart;
    } u;

    LONGLONG QuadPart;
} LARGE_INTEGER;

typedef LARGE_INTEGER *PLARGE_INTEGER;




typedef union _ULARGE_INTEGER {
    struct {
        DWORD LowPart;
        DWORD HighPart;
    } ;
    struct {
        DWORD LowPart;
        DWORD HighPart;
    } u;

    ULONGLONG QuadPart;
} ULARGE_INTEGER;

typedef ULARGE_INTEGER *PULARGE_INTEGER;





typedef LONG_PTR RTL_REFERENCE_COUNT, *PRTL_REFERENCE_COUNT;









typedef struct _LUID {
    DWORD LowPart;
    LONG HighPart;
} LUID, *PLUID;


typedef ULONGLONG  DWORDLONG;
typedef DWORDLONG *PDWORDLONG;





































































































































extern "C" {









unsigned char
__cdecl
_rotl8 (
      unsigned char Value,
      unsigned char Shift
    );

unsigned short
__cdecl
_rotl16 (
      unsigned short Value,
      unsigned char Shift
    );

unsigned char
__cdecl
_rotr8 (
      unsigned char Value,
      unsigned char Shift
    );

unsigned short
__cdecl
_rotr16 (
    );








