# C Coding Style


### Line Width ###

Try to use lines of code between 80 and 120 characters long.

### Indentation ###

We use 4 spaces for indentation.

### Braces ###

1. In GNU style, if either side of an if-else statement has braces, both sides should, to match up indentation:
	
	If (condition):
		{
			Foo;
			Bar;
		}
	else:
		{
			Baz;
		}

2. If the single statement covers multiple lines, e.g. for functions with many arguments, and it is followed by else or else if:

	If (condition):
		{
			a_single_statement_with_many_arguments (some_lengthy_argument,
								         another_lengthy_argument,
								           plus_one);
		}	
	else:
		{
			another_single_statement (arg1, arg2);
		}

3. If the condition is composed of many lines:

	if	(condition1 ||
		(condition2 && condition3) ||
		condition4 ||
		(condition5 && (condition6 || condition7))):
		{
			A_simple_statement
		}

### Conditions ###

Do not check boolean values for equality. By using implicit comparisons, the resulting code can be read more like conversational English. Another rationale is that a ‘true’ value may not be necessarily equal to whatever the TRUE macro uses.
The C language uses the value 0 for many purposes. As a numeric value, the end of a string, a null pointer and the FALSE boolean. To make the code clearer, you should write code that highlights the speci?c way 0 is used. So when reading a comparison, it is possible to know the variable type. For boolean variables, an implicit comparison is appropriate because it’s already a logical expression.

### Functions ###

	The argument list must be broken into a new line for each argument, with the argument names right aligned, taking into account pointers:
	
	Void
	My_function (some_type_t	type,
			Another_type_t	*a_pointer,
			Double_ptr_t		**double_pointer)
	{
		…
	}

The alignment also holds when invoking a function without breaking the line length limit.

### Whitespace ###

When declaring a structure type use newlines to separate logical sections of the structure:

	struct _GtkWrapBoxPrivate
	{
		GtkOrientation        orientation;
		GtkWrapAllocationMode mode;
	
		GtkWrapBoxSpreading   horizontal_spreading;
		GtkWrapBoxSpreading   vertical_spreading;
	}

### The switch Statement ####

A switch should open a block on a new indentation level, and each case should start on the same indentation level as the curly braces, with the case block on a new indentation level:

	switch (condition) {
	case FOO:
		do_foo ();
		break;
	case BAR:
		do_bar;
		break;
	}

### Header Files ###

	return_type	function_name	(type   argument,
                                         type   argument,
                                         type   argument);


The maximum width of each column is given by the longest element in the column:

	void	gtk_type_set_property (GtkType      *type,
                                       const gchar  *value,
                                       GError      **error);
	const gchar *gtk_type_get_property (GtkType      *type);
