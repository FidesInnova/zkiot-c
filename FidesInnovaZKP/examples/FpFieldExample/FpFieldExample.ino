#include <FpField.h>

void setup() {
    Serial.begin(115200);
    delay(500);
    
    // Set the prime number P
    uint64_t prime = 15 * (1UL << 27) + 1;
    Fp::setP(prime);

    Fp a = 123456789;
    Fp b = 987654321;

    Serial.print("a: ");
    a.print(Serial);
    Serial.println();

    Serial.print("b: ");
    b.print(Serial);
    Serial.println();

    Fp c = a + b;
    Serial.print("a + b: ");
    c.print(Serial);
    Serial.println();

    Fp d = a - b;
    Serial.print("a - b: ");
    d.print(Serial);
    Serial.println();

    Fp e = a * b;
    Serial.print("a * b: ");
    e.print(Serial);
    Serial.println();

    Fp f = a / b;
    Serial.print("a / b: ");
    f.print(Serial);
    Serial.println();
}

void loop() {
    // put your main code here, to run repeatedly:
}
