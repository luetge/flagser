#!/usr/bin/env python
from lark import Lark, InlineTransformer, Tree
from functools import reduce
import string

# TODO:
# * Allow the filtration functions to receive parameters

# Why is the Token class not exported by lark?

Str = type(u'')


class Token(str):
    def __new__(cls, type_, value):
        inst = str.__new__(cls, value)
        inst.type = type_
        inst.value = value
        return inst


math_parser = Lark(r"""
    filtrations: filtration+
    filtration: NAME filtration_tag* "=" (expr | cases)
    cases: "cases:" case+
    case: "dim" INTEGER ":" (expr | error) | ELSE ":" (expr | error)
    ELSE: "else"

    expr: term ((ADD | SUBTRACT) term)*
    term: factor ((MULTIPLY | DIVIDE) factor)*
    factor: (base (POWER exponent)?) | (EXP "(" exponent ")")
    base: "(" expr ")" | SIGNED_NUMBER | identifier | function | VAR | array_element
    exponent: "(" expr ")" | SIGNED_NUMBER | identifier | function | VAR | array_element
    error: "error" ESCAPED_STRING

    ADD : "+"
    SUBTRACT : "-"
    MULTIPLY : "*"
    DIVIDE : "/"
    POWER : "^"
    EXP : "exp"

    identifier: DIMENSION | arrays
    arrays: FACE_WEIGHTS | CELL_VERTICES
    global_arrays: VERTEX_WEIGHTS | VERTEX_OUT_DEGREES | VERTEX_IN_DEGREES
    DIMENSION : "dimension"
    FACE_WEIGHTS : "faceWeights"
    CELL_VERTICES : "cellVertices"
    VERTEX_OUT_DEGREES: "vertexOutDegrees"
    VERTEX_IN_DEGREES: "vertexInDegrees"
    VERTEX_WEIGHTS: "vertexWeights"
    filtration_tag: OVERRIDE_VERTICES | OVERRIDE_EDGES
    OVERRIDE_VERTICES: "overrideVertices"
    OVERRIDE_EDGES: "overrideEdges"

    function: max | min | sum | product | map | reduce
    array_functions: map
    array_like: arrays | array_functions
    max: "max(" array_like ("," range)? ")"
    min: "min(" array_like ("," range)? ")"
    sum: "sum(" array_like ("," range)? ("," lambda)? ")"
    product: "product(" array_like ("," range)? ("," lambda)? ")"
    reduce: ("combine" | "reduce") "(" array_like ("," range)? ("," lambda2) ("," expr) ")"
    range: (from "," to)
    map: ("map" | "modifyEach") "(" array_like ("," range)? ("," lambda) ")"
    lambda: VAR "->" expr
    lambda2: VAR VAR "->" expr
    from: expr
    to: expr

    array_element: (array_like | global_arrays) "[" expr "]"

    NAME: /[a-zA-Z_0-9]+/
    VAR: /[a-zA-Z]+/
    INTEGER : /[0-9]+/
    %import common.ESCAPED_STRING
    %import common.SIGNED_NUMBER
    %import common.WS
    %ignore WS

    """, start='filtrations')


def print_error_message(expression):
    return """
            throw std::runtime_error({error_msg});
    """.format(error_msg=expression.children[0])


def generate_code(parse_tree):
    ret = """// WARNING: This code is generated, all changes will be lost after the next compilation.
#ifndef FLAGSER_ALGORITHMS_H
#define FLAGSER_ALGORITHMS_H
#include <string>
#include <vector>

#include "filtration_algorithms.h"
    """
    for filtration in parse_tree.children:
        result = ''
        if filtration.children[-1].data == 'cases':
            for case in filtration.children[-1].children:
                if case.children[0].type == 'ELSE':
                    if hasattr(case.children[1], 'data') and case.children[1].data == 'error':
                        result += print_error_message(
                            case.children[1])
                    else:
                        result += "return {code};".format(code=expression_to_string(
                            case.children[1]))
                else:
                    code = ''
                    if hasattr(case.children[1], 'data') and case.children[1].data == 'error':
                        code = print_error_message(case.children[1])
                    else:
                        code = "return {code};".format(code=expression_to_string(
                            case.children[1]))
                    result += """
                    if (dimension == {dim}) {{
                        {code}
                    }}
                    """.format(dim=case.children[0], code=code)
        else:
            result = "return {code};".format(code=expression_to_string(
                filtration.children[-1]))

        def find_face_filtration(expression):
            if hasattr(expression, 'children'):
                for child in expression.children:
                    if find_face_filtration(child):
                        return True
                return False

            return expression.type == 'FACE_WEIGHTS'
        uses_the_face_filtration = find_face_filtration(filtration)
        face_filtration_code = "\nvirtual inline bool needs_face_filtration() const { return false; }" if not uses_the_face_filtration else ""
        override_vertices_code = ""
        override_edges_code = ""
        for child in filtration.children:
            if hasattr(child, 'data') and child.data == "filtration_tag":
                if child.children[0].type == "OVERRIDE_VERTICES":
                    override_vertices_code = "\nvirtual inline bool overwrite_vertex_filtration() const { return true; }"
                if child.children[0].type == "OVERRIDE_EDGES":
                    override_edges_code = "\nvirtual inline bool overwrite_edge_filtration() const { return true; }"
        ret += """
struct {name}_filtration : public filtration_algorithm_t {{
    virtual inline value_t compute_filtration(unsigned short dimension, const directed_flag_complex_cell_t& cell,
                                    const filtered_directed_graph_t& graph,
                                    const value_t* boundary_filtration) const {{
        {expression}
    }}{face_filtration}{override_vertices}{override_edges}
}};
""".format(
            name=filtration.children[0],
            expression=result,
            face_filtration=face_filtration_code,
            override_vertices=override_vertices_code,
            override_edges=override_edges_code
        )

    # Now add a function that resolves the correct filtration functions
    algorithm_selection_code = "\n".join(list(map(
        lambda expr: "if (algorithm == \"{name}\") return new {name}_filtration();".format(
            name=expr.children[0]),
        parse_tree.children
    )))
    algorithm_list = ", ".join(list(map(
        lambda expr: "\"{name}\"".format(name=expr.children[0]),
        parse_tree.children
    )))
    ret += """
    filtration_algorithm_t* get_custom_filtration_computer(std::string algorithm) {{
{algorithms}

        return nullptr;
    }}
    std::vector<std::string> custom_filtration_computer{{ {algorithm_list} }};
    #endif // FLAGSER_ALGORITHMS_H
    """.format(algorithms=algorithm_selection_code, algorithm_list=algorithm_list)
    return ret


def is_token(expression):
    return not hasattr(expression, 'children')


def token_to_string(expression, bound_variables):
    if expression.type in ['DIMENSION', 'MULTIPLY', 'DIVIDE', 'ADD', 'SUBTRACT']:
        return str(expression)
    if expression.type == 'SIGNED_NUMBER':
        return ((expression + 'f') if '.' in expression else (expression + '.0f'))
    if expression.type == 'VAR':
        # TODO: How can we force identifiers to not be recognised as VARs?
        if expression not in bound_variables and expression != 'dimension':
            raise ValueError("Could not resolve variable " + expression)
        return bound_variables.get(str(expression), str(expression))

    raise ValueError("Invalid token " + str(expression))


def get_array(expression, dive_deeper):
    expr = expression.children[0] if expression.data == 'array_like' else expression
    if expr.data in ['arrays', 'global_arrays']:
        if expr.children[0].type == 'FACE_WEIGHTS':
            return 'boundary_filtration'
        elif expr.children[0].type == 'CELL_VERTICES':
            return 'cell'
        elif expr.children[0].type == 'VERTEX_WEIGHTS':
            return 'graph.vertex_filtration'
    else:
        return dive_deeper(expression)


def range_get_from(range_expression):
    return range_expression.children[0]


def range_get_to(range_expression):
    return range_expression.children[1]


def simple_array_operation(expression, dive_deeper):
    frm = '0' if len(expression.children) == 1 else dive_deeper(
        range_get_from(expression.children[1]))
    to = 'dimension' if len(expression.children) == 1 else dive_deeper(
        range_get_to(expression.children[1]))

    return "{op}({arr}, {frm}, {to})".format(
        op=expression.data,
        arr=get_array(expression.children[0], dive_deeper),
        frm=frm,
        to=to
    )


# The lambda operation always comes last, so check for existence
def is_simple_array_operation(expression):
    return expression.children[-1].data != 'lambda'


def map_(expression, dive_deeper, depth, bound_variables):
    frm = '0' if len(expression.children) == 2 else dive_deeper(
        range_get_from(expression.children[1]))
    to = 'dimension' if len(expression.children) == 2 else dive_deeper(
        range_get_to(expression.children[1]))

    lambda_exp = expression.children[-1]
    counter_name = 'COUNTER_' + (string.ascii_uppercase[depth])

    # Compute the formula for the lambda function
    new_bound_variables = bound_variables.copy()
    new_bound_variables.update(
        {str(lambda_exp.children[0]): get_array(expression.children[0], dive_deeper) + '[' + counter_name + ']'})
    value = dive_deeper(lambda_exp, depth + 1, new_bound_variables)

    # Execute code as IIFE
    return """
        [&](){{
            std::vector<value_t> EXPRESSION_{depth};
            for (int {counter} = (int){frm}; {counter} <= (int){to}; {counter}++) {{
                EXPRESSION_{depth}.push_back({value});
            }}
            return EXPRESSION_{depth};
        }}()""".format(
        counter=counter_name,
        frm=frm,
        to=to,
        value=value,
        depth=depth
    )


def reduce_(expression, dive_deeper, depth, bound_variables, operation='='):
    frm = '0' if len(expression.children) == 3 else dive_deeper(
        range_get_from(expression.children[1]))
    to = 'dimension' if len(expression.children) == 3 else dive_deeper(
        range_get_to(expression.children[1]))

    lambda_exp = expression.children[-2]
    counter_name = 'COUNTER_' + (string.ascii_uppercase[depth])

    # Compute the formula for the lambda function
    new_bound_variables = bound_variables.copy()
    new_bound_variables.update({
        str(lambda_exp.children[0]): 'EXPRESSION_{depth}'.format(depth=depth),
        str(lambda_exp.children[1]): get_array(expression.children[0], dive_deeper) + '[' + counter_name + ']',
    })
    formula = dive_deeper(lambda_exp, depth + 1, new_bound_variables)

    # Execute code as IIFE
    return """
        [&](){{
            value_t EXPRESSION_{depth} = {initial_value};
            for (int {counter} = (int){frm}; {counter} <= (int){to}; {counter}++) {{
                EXPRESSION_{depth} {operation} {formula};
            }}
            return EXPRESSION_{depth};
        }}()""".format(
        initial_value=dive_deeper(
            expression.children[-1], depth, bound_variables),
        counter=counter_name,
        frm=frm,
        to=to,
        formula=formula,
        operation=operation,
        depth=depth
    )


def expression_to_string(expression, depth=0, bound_variables={}):
    def dive_deeper(exp, d=depth, bV=bound_variables.copy()):
        return expression_to_string(exp, d, bV)

    def continue_with_children(exp, d=depth, bV=bound_variables.copy()):
        if len(expression.children) > 1:
            return reduce(lambda carry, exp: carry + dive_deeper(exp), expression.children, '')
        return dive_deeper(expression.children[0])

    # If we are at leaves, handle the corresponding expressions
    if is_token(expression):
        return token_to_string(expression, bound_variables)

    # In this case, we have a tree element consisting of data (the identifier) and children
    if expression.data in ['max', 'min', 'sum', 'product']:
        if is_simple_array_operation(expression):
            return simple_array_operation(expression, dive_deeper)

        # Here we can only have sum or product with lambda
        operation = '+=' if expression.data == 'sum' else '*='
        base_value = '0' if expression.data == 'sum' else '1'
        modified_tree = Tree('reduce', expression.children[0:-1] + [Tree('lambda', [
            Token('VAR', '_c_a_r_r_y')] + expression.children[-1].children), Token('SIGNED_NUMBER', base_value)])
        return reduce_(modified_tree, dive_deeper, depth, bound_variables, operation)

    if expression.data == 'map':
        return map_(expression, dive_deeper, depth, bound_variables)

    if expression.data == 'reduce':
        return reduce_(expression, dive_deeper, depth, bound_variables)

    if expression.data == 'factor':
        if hasattr(expression.children[0], 'type') and expression.children[0].type == 'EXP':
            return 'exp({expr})'.format(expr=dive_deeper(expression.children[1]))
        if hasattr(expression.children[0], 'data') and expression.children[0].data == 'base' and len(expression.children) > 1:
            return 'pow({base}, {exp})'.format(base=dive_deeper(expression.children[0]), exp=dive_deeper(expression.children[2]))

    if expression.data == 'array_element':
        if hasattr(expression.children[0], 'data') and expression.children[0].data == 'global_arrays':
            if expression.children[0].children[0].type == 'VERTEX_OUT_DEGREES':
                return "graph.outdegrees[(int)({idx})]".format(idx=dive_deeper(expression.children[1]))
            if expression.children[0].children[0].type == 'VERTEX_IN_DEGREES':
                return "graph.indegrees[(int)({idx})]".format(idx=dive_deeper(expression.children[1]))

        return "{arr}[(int)({idx})]".format(
            arr=get_array(expression.children[0], dive_deeper),
            idx=dive_deeper(expression.children[1])
        )

    if expression.data in ['lambda', 'lambda2']:
        # Just evaluate the function itself without the variable declaration
        return dive_deeper(expression.children[-1])

    if expression.data == 'expr':
        return "({expr})".format(expr=continue_with_children(expression))

    return continue_with_children(expression)


if __name__ == "__main__":
    import os
    dir_path = os.path.dirname(os.path.realpath(__file__))
    file = open(os.path.join(dir_path, '../algorithms.math'), 'r')
    code = generate_code(math_parser.parse(file.read()))
    file.close()
    of = open(os.path.join(dir_path, '../include/algorithms.h'), 'w')
    of.write(code)
    of.close()
