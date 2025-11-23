/**
 * ErrorBoundary - React 错误边界组件
 * 捕获子组件错误并显示友好的错误界面
 */
import { Component } from 'react'
import type { ErrorInfo, ReactNode } from 'react'
import { Result, Button, Typography, Card } from 'antd'
import { ExclamationCircleOutlined, ReloadOutlined, HomeOutlined } from '@ant-design/icons'

const { Paragraph, Text } = Typography

interface ErrorBoundaryProps {
  children: ReactNode
  fallback?: ReactNode
  onReset?: () => void
}

interface ErrorBoundaryState {
  hasError: boolean
  error: Error | null
  errorInfo: ErrorInfo | null
}

export class ErrorBoundary extends Component<ErrorBoundaryProps, ErrorBoundaryState> {
  constructor(props: ErrorBoundaryProps) {
    super(props)
    this.state = {
      hasError: false,
      error: null,
      errorInfo: null,
    }
  }

  static getDerivedStateFromError(error: Error): Partial<ErrorBoundaryState> {
    return { hasError: true, error }
  }

  componentDidCatch(error: Error, errorInfo: ErrorInfo) {
    // 记录错误到日志服务
    console.error('ErrorBoundary caught an error:', error, errorInfo)

    this.setState({
      error,
      errorInfo,
    })
  }

  handleReset = () => {
    this.setState({
      hasError: false,
      error: null,
      errorInfo: null,
    })
    this.props.onReset?.()
  }

  handleGoHome = () => {
    window.location.href = '/'
  }

  render() {
    if (this.state.hasError) {
      // 使用自定义 fallback
      if (this.props.fallback) {
        return this.props.fallback
      }

      // 默认错误界面
      return (
        <div
          style={{
            minHeight: '100vh',
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'center',
            padding: '24px',
            background: 'var(--bg-secondary)',
          }}
        >
          <Card style={{ maxWidth: 600, width: '100%' }}>
            <Result
              status="error"
              icon={<ExclamationCircleOutlined />}
              title="Oops! Something went wrong"
              subTitle="An unexpected error occurred. Please try again or contact support if the problem persists."
              extra={[
                <Button type="primary" icon={<ReloadOutlined />} onClick={this.handleReset} key="reset">
                  Try Again
                </Button>,
                <Button icon={<HomeOutlined />} onClick={this.handleGoHome} key="home">
                  Go Home
                </Button>,
              ]}
            >
              {process.env.NODE_ENV === 'development' && this.state.error && (
                <div style={{ marginTop: 24 }}>
                  <Paragraph>
                    <Text strong style={{ fontSize: 16 }}>
                      Error Details:
                    </Text>
                  </Paragraph>
                  <Card
                    size="small"
                    style={{
                      background: 'var(--bg-tertiary)',
                      marginTop: 12,
                    }}
                  >
                    <pre
                      style={{
                        margin: 0,
                        fontSize: 12,
                        fontFamily: 'var(--font-family-mono)',
                        color: 'var(--color-error)',
                        overflow: 'auto',
                        maxHeight: 200,
                      }}
                    >
                      {this.state.error.toString()}
                      {'\n\n'}
                      {this.state.errorInfo?.componentStack}
                    </pre>
                  </Card>
                </div>
              )}
            </Result>
          </Card>
        </div>
      )
    }

    return this.props.children
  }
}

export default ErrorBoundary
