#ifdef CRAYFISH
#  define cforcevector cdir$ ivdep
#else
#    ifdef IBM3090
#      define cforcevector CIGNR IGNORE RECRDEPS 
@PROCESS DIRECTIVE ('IGNR')
#    else
#        ifdef UXP
#          define cforcevector *vocl loop,novrec
#          define cforcescalar *vocl loop,scalar
#          define cforcenoeval !OCL NOEVAL
#          define cforceeval !OCL EVAL
#        else
#          ifdef NEC_SX
#            define cforcevector *vdir nodep
#          else
#            define cforcevector c$dir no_recurrence
#          endif
#        endif
#    endif
#endif
