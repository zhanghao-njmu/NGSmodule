# Multi-stage build for production

# Stage 1: Build
FROM node:18-alpine as build

WORKDIR /app

# Copy package files
COPY package*.json ./

# Install dependencies
RUN npm ci

# Copy source code
COPY . .

# Build arguments
ARG VITE_API_URL
ARG VITE_WS_URL

# Set environment variables for build
ENV VITE_API_URL=${VITE_API_URL}
ENV VITE_WS_URL=${VITE_WS_URL}

# Build application
RUN npm run build

# Stage 2: Production
FROM nginx:alpine

# Copy built files from build stage
COPY --from=build /app/dist /usr/share/nginx/html

# Copy custom nginx config
COPY nginx.conf /etc/nginx/conf.d/default.conf

# Create non-root user
RUN adduser -D -u 1000 nginx-user && \
    chown -R nginx-user:nginx-user /usr/share/nginx/html && \
    chown -R nginx-user:nginx-user /var/cache/nginx && \
    chown -R nginx-user:nginx-user /var/log/nginx && \
    touch /var/run/nginx.pid && \
    chown -R nginx-user:nginx-user /var/run/nginx.pid

# Switch to non-root user
USER nginx-user

# Expose port
EXPOSE 80

# Health check
# Use 127.0.0.1 instead of localhost: in this image, BusyBox wget resolves
# `localhost` to IPv6 `::1`, but nginx only listens on IPv4 0.0.0.0:80,
# so the resolver-IPv6-first behaviour causes a permanent ECONNREFUSED.
# /health is defined in nginx.conf and returns 200 "healthy".
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD wget --quiet --tries=1 --spider http://127.0.0.1/health || exit 1

# Start nginx
CMD ["nginx", "-g", "daemon off;"]
