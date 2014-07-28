#### Multiple Observers

Multiple observers can be used within a simulation. By default the events detected by every observer will written to every output file. To set an individual output file for each observer use flag and condition as follows:

```python
obs1 = SmallObserverSphere(pos1, radius, flag='seen_by_obs1')
obs2 = SmallObserverSphere(pos2, radius, flag='seen_by_obs2')
out1 = ConditionalOutput('events1.txt', condition='seen_by_obs1')
out2 = ConditionalOutput('events2.txt', condition='seen_by_obs2')
```