# Rate Limiting Documentation

## Overview

NGSmodule API implements rate limiting to protect against abuse, DDoS attacks, and ensure fair resource allocation among users.

## Implementation

**Library**: `slowapi` (FastAPI port of Flask-Limiter)
**Storage**: In-memory (default) or Redis (production recommended)
**Strategy**: IP-based for anonymous requests, User-based for authenticated requests

## Rate Limit Configurations

### Global Default
- **Limit**: 100 requests per minute
- **Applies to**: All endpoints by default

### Authentication Endpoints

| Endpoint | Limit | Reason |
|----------|-------|--------|
| `POST /auth/login` | 5/min | Prevent brute-force attacks |
| `POST /auth/register` | 3/min | Prevent spam registrations |

### File Operations

| Endpoint | Limit | Reason |
|----------|-------|--------|
| `POST /files/upload` | 10/min | Prevent storage abuse |
| `GET /files/{id}/download` | 30/min | Balance server load |
| `DELETE /files/{id}` | 20/min | Prevent accidental mass deletion |

### Batch Operations

| Endpoint | Limit | Reason |
|----------|-------|--------|
| `POST /samples/import-csv` | 5/min | Resource-intensive operation |

## Response Format

When rate limit is exceeded, the API returns:

```json
{
  "error": "Rate limit exceeded",
  "detail": "Too many requests. Please slow down and try again later.",
  "code": "RATE_LIMIT_EXCEEDED",
  "retry_after": 60
}
```

**HTTP Status Code**: `429 Too Many Requests`
**Headers**: `Retry-After: <seconds>`

## Configuration

### Using In-Memory Storage (Development)

Default configuration in `app/core/rate_limit.py`:

```python
limiter = Limiter(
    key_func=get_remote_address,
    default_limits=["100/minute"],
    storage_uri="memory://"
)
```

### Using Redis Storage (Production)

Update `.env` file:

```bash
REDIS_URL=redis://localhost:6379/1
```

Then in `app/main.py` startup event:

```python
from app.core.rate_limit import configure_redis_storage

@app.on_event("startup")
async def startup_event():
    init_db()
    configure_redis_storage(settings.REDIS_URL)
```

## Customizing Rate Limits

### Add Rate Limit to New Endpoint

1. Import dependencies:
```python
from fastapi import Request
from app.core.rate_limit import limiter, RateLimits
```

2. Add decorator and Request parameter:
```python
@router.post("/endpoint")
@limiter.limit(RateLimits.CREATE)  # or custom limit like "50/minute"
async def my_endpoint(
    request: Request,  # Required for rate limiting
    # ... other parameters
):
    pass
```

### Define Custom Rate Limit

Add to `RateLimits` class in `app/core/rate_limit.py`:

```python
class RateLimits:
    MY_CUSTOM_LIMIT = "25/minute"
```

## Monitoring

Rate limit events are logged with:
- Client IP address
- Requested endpoint
- Timestamp

Check logs for patterns:
```bash
grep "Rate limit exceeded" logs/app.log
```

## Testing

Simulate rate limit exceeded:

```bash
# Make 6 rapid login attempts (limit is 5/min)
for i in {1..6}; do
  curl -X POST http://localhost:8000/api/v1/auth/login \
    -d "username=test&password=test"
done
```

Expected 6th response: `429 Too Many Requests`

## Whitelist IPs (Future Enhancement)

To exempt specific IPs from rate limiting, modify `get_rate_limit_key()` in `rate_limit.py`:

```python
WHITELISTED_IPS = ["192.168.1.100", "10.0.0.1"]

def get_rate_limit_key(request: Request) -> str:
    ip = get_remote_address(request)
    if ip in WHITELISTED_IPS:
        return "whitelist"  # Shared key for unlimited access
    # ... existing logic
```

## Security Best Practices

1. **Always use Redis in production** - In-memory storage is lost on restart
2. **Monitor rate limit violations** - Frequent violations may indicate attack
3. **Adjust limits based on usage patterns** - Review metrics regularly
4. **Use HTTPS** - Prevent header spoofing
5. **Consider user-based limits** - More accurate for authenticated endpoints

## Troubleshooting

### Rate Limit Not Working

1. Verify `app.state.limiter` is set in `main.py`
2. Check `Request` parameter is first in endpoint
3. Ensure decorator is directly above function

### False Positives

If legitimate users hit limits:
1. Review limit values in `RateLimits` class
2. Check if multiple users share same IP (NAT/proxy)
3. Consider switching to user-based limits for authenticated endpoints

## References

- [slowapi Documentation](https://slowapi.readthedocs.io/)
- [OWASP Rate Limiting Guidelines](https://cheatsheetseries.owasp.org/cheatsheets/Denial_of_Service_Cheat_Sheet.html)
