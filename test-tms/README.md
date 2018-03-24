Tests
=======
We use [Catch2](https://github.com/catchorg/Catch2) as testing framework. [I](https://github.com/catchorg/Catch2) choose it because it's very, very simple. Just one header file.

## Writing Tests

You can start with simple example.

This is a function to be tested:

```c++
unsigned int Factorial( unsigned int number ) {
    return number <= 1 ? number : Factorial(number-1)*number;
}
```

This is a compited test:

```c++
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"

unsigned int Factorial( unsigned int number ) {
    return number <= 1 ? number : Factorial(number-1)*number;
}

TEST_CASE( "Factorials are computed", "[factorial]" ) {
    REQUIRE( Factorial(1) == 1 );
    REQUIRE( Factorial(2) == 2 );
    REQUIRE( Factorial(3) == 6 );
    REQUIRE( Factorial(10) == 3628800 );
}
```

It's all! Very simple!

You can [see more documentation and examples](https://github.com/catchorg/Catch2/tree/master/docs#reference).