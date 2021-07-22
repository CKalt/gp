pub mod parsers {
    use nom::{
        bytes::complete::take_while,
        character::complete::{alphanumeric1, char},
        IResult,
        branch::alt,
        sequence::{preceded, delimited, pair},
        multi::many1,
    };

    #[derive(Debug)]
    pub enum ParsedSNode<'a> {
        Term(&'a str),
        Func(&'a str, Vec<ParsedSNode<'a>>),
    }
    use ParsedSNode::*;

    pub type ParsedArgs<'a> = Vec<ParsedSNode<'a>>;

    pub fn parse_sexpr(input: &str) -> IResult<&str, ParsedSNode> {
        alt(
            (
                parse_term_node,
                parse_func_node,
            )
        )(input)
    }

    pub fn parse_term_node(input: &str) -> IResult<&str, ParsedSNode> {
        preceded(
            skip_spaces,
            alphanumeric1,
        )(input).map(|(input, output)| (input, Term(output)))
    }

    pub fn parse_func_node(input: &str) -> IResult<&str, ParsedSNode> {
        preceded(
            skip_spaces,
            delimited(
                char('('),
                parse_func_call, 
                preceded(
                    skip_spaces,
                    char(')')
                )
            )
        )(input)
    }

    pub fn parse_func_call(input: &str) -> IResult<&str, ParsedSNode> {
        pair(parse_func_name, parse_func_args)(input).map(|(input, output)|
            (input, Func(output.0, output.1)))
    }

    pub fn parse_func_name(input: &str) -> IResult<&str, &str> {
        preceded(
            skip_spaces,
            alphanumeric1,
        )(input)
    }

    pub fn parse_func_args(input: &str) -> IResult<&str, ParsedArgs> {
        many1(
            alt(
                (
                    parse_term_node,
                    parse_func_node,
                )
            )
        )(input)
    }

    pub fn skip_spaces(input: &str) -> IResult<&str, &str> {
        let chars = " \t\r\n";
        take_while(move |ch| chars.contains(ch))(input)
    }

    #[cfg(test)]
    mod tests {
            use super::*;

            #[test]
            fn test_parse_sexpr() {
                assert_eq!(
                    format!("{:?}", parse_sexpr("(ADD X Y)")), 
                            r#"Ok(("", Func("ADD", [Term("X"), Term("Y")])))"#);
                assert_eq!(
                    format!("{:?}", parse_sexpr("( ADD  X    Y )"
                            )), 
                            r#"Ok(("", Func("ADD", [Term("X"), Term("Y")])))"#);

                assert_eq!(
                    format!("{:?}", parse_sexpr(
                        "( ADD    X (DIV (IF 3 1 3)  2) ( MULT 1 4 1  )  )"
                            )), 
                            r#"Ok(("", Func("ADD", [Term("X"), Func("DIV", [Func("IF", [Term("3"), Term("1"), Term("3")]), Term("2")]), Func("MULT", [Term("1"), Term("4"), Term("1")])])))"#);

                assert_eq!(
                    format!("{:?}", parse_sexpr("( ADD ()   X Y  )")),
                            r#"Err(Error(Error { input: ")   X Y  )", code: AlphaNumeric }))"#);

            }

            #[test]
            fn test_parse_term_node() {
                assert_eq!(
                    format!("{:?}", parse_term_node("X")), 
                            r#"Ok(("", Term("X")))"#);
                assert_eq!(
                    format!("{:?}", parse_term_node("X Y Z")), 
                            r#"Ok((" Y Z", Term("X")))"#);
                assert_eq!(
                    format!("{:?}", parse_term_node("  X123YZ Y92 Z29 ")),
                            r#"Ok((" Y92 Z29 ", Term("X123YZ")))"#);
                assert_eq!(
                    format!("{:?}", parse_term_node(")  X123YZ Y92 Z29 ")),
                            r#"Err(Error(Error { input: ")  X123YZ Y92 Z29 ", code: AlphaNumeric }))"#);
            }

            #[test]
            fn test_parse_func_node() {
                assert_eq!(
                    format!("{:?}", parse_func_node(" ( XYZ 1  2   (AD 3) ) ")),
                            r#"Ok((" ", Func("XYZ", [Term("1"), Term("2"), Func("AD", [Term("3")])])))"#);

                assert_eq!(
                    format!("{:?}", parse_func_node(" ( XYZ 1  2   (AD 3 ) ")),
                            r#"Err(Error(Error { input: "", code: Char }))"#);
            }

            #[test]
            fn test_parse_func_call() {
                assert_eq!(
                    format!("{:?}", parse_func_call("  ABCD123 X1 (Y 1) Z  ")),
                            r#"Ok(("  ", Func("ABCD123", [Term("X1"), Func("Y", [Term("1")]), Term("Z")])))"#);

                assert_eq!(
                    format!("{:?}", parse_func_call("  (ABCD123 X1 (Y 1) Z)  ")),
                            r#"Err(Error(Error { input: "(ABCD123 X1 (Y 1) Z)  ", code: AlphaNumeric }))"#);

            }

            #[test]
            fn test_parse_func_name() {
                assert_eq!(
                    format!("{:?}", parse_func_name("    ABCD123 ")),
                            r#"Ok((" ", "ABCD123"))"#);
                assert_eq!(
                    format!("{:?}", parse_func_name(" (   ABCD123 )")),
                            r#"Err(Error(Error { input: "(   ABCD123 )", code: AlphaNumeric }))"#);
            }

            #[test]
            fn test_parse_func_args() {
                assert_eq!(
                    format!("{:?}", parse_func_args("  X1 (Y 1) Z  ")),
                            r#"Ok(("  ", [Term("X1"), Func("Y", [Term("1")]), Term("Z")]))"#);
                assert_eq!(
                    format!("{:?}", parse_func_args("  (X1 (Y 1) Z   ")),
                            r#"Err(Error(Error { input: "", code: Char }))"#);

                assert_eq!(
                    format!("{:?}", parse_func_args("  ((X1 Y) X) (Y 1) Z   ")),
                            r#"Err(Error(Error { input: "(X1 Y) X) (Y 1) Z   ", code: AlphaNumeric }))"#);
              }

            #[test]
            fn test_skip_spaces() {
                assert_eq!(
                    format!("{:?}", skip_spaces("  \t \r  \n\n AABC DEF ")),
                            r#"Ok(("AABC DEF ", "  \t \r  \n\n "))"#);
            }
    }
}
