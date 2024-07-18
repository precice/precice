These tests couple 4 participants:

```
EA - IA - IB - EB
```

* EA and IA are explicitly coupled.
* EB and IB are explicitly coupled.
* IA and IB are implicitly coupled.

The tests in the suites `implicit-first` and `implicit-second` run the implicit scheme before the explicit schemes and vice versa.

The tests write data with the format `timewindow.iteration` or `0.1` for initial data.
They sample at the end of the time window and check if the explected sample is present.
