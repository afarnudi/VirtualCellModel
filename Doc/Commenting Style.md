# Commenting Style

## Description Before Members ##

### Detailed Description ###

For each entity there are two (in some cases three) kinds of descriptions, which together form the description for that entity. Brief and detailed descriptions both are optional. For methods and functions there is a third type of descriptions, the so called in body description, which consists of the concatenation of all the comment blocks found within the body od the method or function.
Description Before Members
Detailed Description
Use this method before entity:

        ///
        /// ... text ...
        /// 

### Brief Description ###

Use the \brief command with above comment block. This comment ends at the end of a paragraph, so the detailed description follows after an empty line.
Here is an example:
/*! \brief Brief description.
*         Brief description continued.
*
*  Detailed description starts here.
*/

## Description After Members ##

## Detailed Description ##

	Int var;    ///< detailed description after member.
		    ///<

## Brief description ##

	Int var;		///< brief description after the member

Documenting members manually
For each entity we have a structural command, here is the list of structural commands:

\ class to document a class\*
\ struct to document a C-struct
\ union to document a union
\ enum to document an enumeration type
\ fn to document a function
\ var to document a variable or typedef of enum value
\ def to document a #define
\ typedef to document a type defenition
\ file to document a file
\ namespace to document a namespace
\ package to document a Java package
\ interface to document an IDL interface

Here is an example:

	///	\ class test
	/// 	\brief a class test
	///	
	///	a more detailed class descriptions.
	///

Attention

To document global objects (functions, typedefs, enum, macros, etc), you must document the file in which they are defined. In other words, there must at least be a
	
	/// \ file
	/// 
	///
