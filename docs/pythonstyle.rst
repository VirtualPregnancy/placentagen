============
Python Style
============

In general, one should try to follow the PEP8 style for python programming. Details of the style can be found
`at this link <https://legacy.python.org/dev/peps/pep-0008/>`_.

Some IDEs (integrated development environments) are able to check your code for violation of PEP8 conventions,
and are worth investigating.

The code is split into modules which contain collections of functions that do related things. Each module can be
automatically documented to some extent if you follow conventions in inserting `doc strings` into the code.

In general, at a minimum any module or function should contain the following:

.. code-block:: python

    #!/usr/bin/env python
    import whatever.libraries.you.need

    """
    .. module:: module_name
        :synopsis: One sentence synopsis

    :synopsis:A longer synopsis that could appear on the home page
        for that module in documentation.

    """


    def function_name(your, inputs):
        """ What this function does

        Inputs:
        - input name: A description of the input
        - name 2: another description

        Returns:
            - output name: a description of the output
                - can do sub points if there are complicated arrays


        A way you might want to use me is:

        >>> include a simple example of how to call the function

        Tell us what the function will do in this example case
        """
        some_fun_stuff = 2.0 * inputs / your

        return some_fun_stuff


