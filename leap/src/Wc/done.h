#ifndef _Wc_done_h_
#define _Wc_done_h_

/* -- converter done macro for `old-style' Xrm converters
*******************************************************************************
    One of these two macros are invoked when a resource converter has
    completed the conversion - `done' for old-style converters, and
    `new_done' for new-style converters.  The first is taken directly from
    the Xt Reference Manual which came with Release 4 of the X11 Window
    System from MIT.
*/
#ifndef done

#define done( type, value ) 			\
{						\
    if ( toVal->addr != NULL )			\
    {						\
	if ( toVal->size < sizeof( type ) )	\
	{					\
	    toVal->size = sizeof( type );	\
	    return;				\
	}					\
	*(type*)(toVal->addr) = (value);	\
    }						\
    else					\
    {						\
	static type static_val;			\
	static_val = (value);			\
	toVal->addr = (caddr_t)&static_val;	\
    }						\
    toVal->size = sizeof(type);			\
    return;					\
}

#endif /* done */


/* -- converter new_done macro for `new-style' Xrm converters
*******************************************************************************
    This macro is taken from page 789 out of the extremely useful book
    by Paul Asente and Ralph Swick titled "X Window System Toolkit" and
    published by Digital Press, ISBN 1-55558-051-3.
*/
#ifndef new_done

#define new_done( type, value )			\
{						\
    if ( toVal->addr != NULL )			\
    {						\
	if ( toVal->size < sizeof( type ) )	\
	{					\
	    toVal->size = sizeof( type );	\
	    return False;			\
	}					\
	*(type*)(toVal->addr) = (value);	\
    }						\
    else					\
    {						\
	static type static_val;			\
	static_val = (value);			\
	toVal->addr = (caddr_t)&static_val;	\
    }						\
    toVal->size = sizeof(type);			\
    return True;				\
}

#endif /* new_done */

#endif /* _Wc_done_h_ */
